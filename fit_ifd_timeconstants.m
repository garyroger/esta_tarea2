function [Tprime, Tpp, R] = fit_ifd_timeconstants(t, y, t0, Tw)
% FIT_IFD_TIMECONSTANTS  estima T'_{do} y T''_{do} de i_fd(t)
% usa ajuste bi-exponencial y metodo grafico 0.368*h (ambos).
%
% [Tprime, Tpp, R] = fit_ifd_timeconstants(t, y, t0, Tw)
%   t,y : vectores de tiempo y i_fd (misma longitud)
%   t0  : instante del evento (seg)
%   Tw  : (opcional) ventana de analisis desde t0 (seg). default: cubre todo.
%
% Salidas:
%   Tprime = T'_{do} (transitoria, s)
%   Tpp    = T''_{do} (subtransitoria, s)
%   R      = struct con parametros, ajustes y resultados 0.368*h

t = t(:); y = y(:);
assert(numel(t)==numel(y), 't e y deben tener misma longitud');

% ventana de analisis
if nargin<4 || isempty(Tw), Tw = t(end)-t0; end
idx = (t>=t0) & (t<=t0+Tw);
tt = t(idx); yy = y(idx);

% estimar asintotico (promedio final del 20% ultimo de la ventana)
tailStart = tt(1) + 0.8*(tt(end)-tt(1));
Iinf0 = mean(yy(tt>=tailStart));

% semillas
A1_0 = yy(find(tt>=tailStart,1)) - Iinf0;   % cola
A1_0 = max(min(A1_0,  2), -2);              % limitar
A2_0 = (yy(1) - Iinf0) - A1_0;              % resto inicial
T1_0 = 3.0;                                 % s (transitorio)
T2_0 = 0.03;                                % s (subtransitorio)
theta0 = [Iinf0, A1_0, T1_0, A2_0, T2_0];

% modelo
f = @(th, x) th(1) + th(2).*exp(-(x-t0)./max(th(3),1e-3)) + ...
                    th(4).*exp(-(x-t0)./max(th(5),1e-4));

% ajuste (lsqcurvefit si existe; sino fminsearch)
use_lsq = exist('lsqcurvefit','file')==2;
if use_lsq
    lb = [Iinf0-0.2, -5, 0.5, -5, 0.005];
    ub = [Iinf0+0.2,  5, 10,  5, 0.2  ];
    opts = optimoptions('lsqcurvefit','Display','off');
    theta = lsqcurvefit(@(th,x) f(th,x), theta0, tt, yy, lb, ub, opts);
else
    % fminsearch en espacio log para tiempos (para mantener >0)
    g = @(p) f([p(1), p(2), exp(p(3)), p(4), exp(p(5))], tt) - yy;
    p0 = [theta0(1) theta0(2) log(theta0(3)) theta0(4) log(theta0(5))];
    p  = fminsearch(@(pp) norm(g(pp)), p0, optimset('Display','off'));
    theta = [p(1), p(2), exp(p(3)), p(4), exp(p(5))];
end

Iinf   = theta(1);
A1     = theta(2);  Tprime = theta(3);
A2     = theta(4);  Tpp    = theta(5);

% --- metodo 0.368*h ---
% subtransitorio: desde t0, comparando a Iinf con la parte temprana
h1 = abs(yy(1) - Iinf);
y0368_pp = Iinf + sign(yy(1)-Iinf)*0.368*h1;
t0368_pp = interp_time(tt, yy, y0368_pp, [t0, t0+min(0.2, 4*Tpp)]);

% transitorio: inicia tras morir el rapido (t_tail = t0 + 4*Tpp)
t_tail   = t0 + 4*Tpp;
% valor al inicio de la cola:
y_tail   = interp1(tt, yy, t_tail, 'linear', 'extrap');
h2       = abs(y_tail - Iinf);
y0368_p  = Iinf + sign(y_tail-Iinf)*0.368*h2;
t0368_p  = interp_time(tt, yy, y0368_p, [t_tail, t_tail+max(1.0, 2*Tprime)]);

Tpp_0368 = max(t0368_pp - t0, 0);
Tp_0368  = max(t0368_p  - t_tail, 0);

% --- graficas ---
figure; clf
subplot(1,2,1)
plot(t,y,'k','LineWidth',1.2); hold on; grid on
plot(tt, f(theta,tt),'r--','LineWidth',1.2)
xline(t0,':'); yline(Iinf,':','I_\infty');
title(sprintf('ajuste bi-exponencial: T''=%.4fs, T''''=%.4fs', Tprime, Tpp))
xlabel('t [s]'); ylabel('i_{fd} [pu]')
legend('i_{fd}','ajuste','Location','best')

subplot(1,2,2)
semilogy(tt, abs(yy-Iinf), 'k','LineWidth',1.2); hold on; grid on
% rectas teoricas
semilogy(tt, abs(A1)*exp(-(tt-t0)/Tprime), 'r--', 'LineWidth',1.1)
semilogy(tt, abs(A2)*exp(-(tt-t0)/Tpp),    'b--', 'LineWidth',1.1)
xline(t_tail,':','t_{tail}=t_0+4T'''''); 
title('semilog: |i_{fd}-I_\infty| y componentes')
xlabel('t [s]'); ylabel('|i_{fd}-I_\infty|')

% marcadores 0.368*h
subplot(1,2,1); hold on
plot([t0 t0368_pp], [y0368_pp y0368_pp], 'b:','LineWidth',1.1)
plot(t0368_pp, y0368_pp, 'bo','MarkerFaceColor','b')
text(t0368_pp, y0368_pp, sprintf('  T''''=%.3f s (0.368h)',Tpp_0368))

plot([t_tail t0368_p], [y0368_p y0368_p], 'm:','LineWidth',1.1)
plot(t0368_p, y0368_p, 'mo','MarkerFaceColor','m')
text(t0368_p, y0368_p, sprintf('  T''=%.3f s (0.368h)',Tp_0368))

% --- reporte ---
fprintf('=== i_fd: estimacion de constantes de tiempo ===\n');
fprintf('Iinf       = %.5f pu\n', Iinf);
fprintf('A1 (trans) = %.5f,  T''do  (ajuste) = %.4f s\n', A1, Tprime);
fprintf('A2 (subtr) = %.5f,  T''''do (ajuste) = %.4f s\n', A2, Tpp);
fprintf('T''do  (0.368h) = %.4f s  |  T''''do (0.368h) = %.4f s\n', Tp_0368, Tpp_0368);

% salida extra
R = struct('Iinf',Iinf,'A1',A1,'A2',A2,'Tprime_fit',Tprime,'Tpp_fit',Tpp, ...
           'Tprime_0368',Tp_0368,'Tpp_0368',Tpp_0368, ...
           't0',t0,'t_tail',t_tail,'theta',theta);

end

% --- helper: encontrar t donde y cruza yref en una ventana [ta,tb]
function tc = interp_time(t, y, yref, win)
if nargin<4, win = [t(1) t(end)]; end
mask = (t>=win(1) & t<=win(2));
t2 = t(mask); y2 = y(mask);
% buscamos cruce por diferencia minima
[~,k] = min(abs(y2 - yref));
k = max(2,min(k,numel(t2)));
% interpola lineal entre (k-1,k)
tc = interp1(y2(k-1:k), t2(k-1:k), yref, 'linear', 'extrap');
end
