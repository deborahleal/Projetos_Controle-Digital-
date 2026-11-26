clear; clc; close all;

% Questão 2 - Controle Digital
Ts = 0.1;
Ad = [0.9142 0.1910; 0.0955 1.0097];
Bd = [0.0955; 0.0048];
C  = [0 1];
D  = 0;
n  = size(Ad,1);
sysd = ss(Ad, Bd, C, D, Ts);

%% (a) Estabilidade
polos_ma = eig(Ad)
modulos  = abs(polos_ma)
% polo em z = 1,1052, fora do círculo unitário -> sistema instável

%% (b) Controlabilidade
Co = ctrb(Ad, Bd)
det_Co  = det(Co)
posto_Co = rank(Co)

%% (c) Observabilidade
Ob = obsv(Ad, C)
det_Ob  = det(Ob)
posto_Ob = rank(Ob)

%% (d) Diagrama de blocos
% feito no draw.io (u[k] = N r[k] - K x[k])

%% (e) Projeto de K por alocação de polos
xi = 0.7;
wn = 3;
wd = wn*sqrt(1 - xi^2);

s_des = [-xi*wn + 1j*wd; -xi*wn - 1j*wd];
z_des = exp(s_des*Ts)
alpha = poly(z_des)

K   = acker(Ad, Bd, z_des)
Acl = Ad - Bd*K;
polos_mf = eig(Acl)

N = 1/(C/(eye(n) - Acl)*Bd)

%% (f) Função de transferência
Gz = minreal(tf(sysd))
syscl = ss(Acl, Bd*N, C, D, Ts);
Tz = minreal(tf(syscl))

%% (g) Resposta ao degrau
t = 0:Ts:5;
r = ones(size(t));
x = zeros(n, length(t));
u = zeros(size(t));
y = zeros(size(t));

for k = 1:length(t)
    u(k) = N*r(k) - K*x(:,k);
    y(k) = C*x(:,k);
    if k < length(t)
        x(:,k+1) = Ad*x(:,k) + Bd*u(k);
    end
end

y_ss = dcgain(syscl);
e_ss = 1 - y_ss

[yp, ip] = max(y);
Mp = 100*(yp - y_ss)/y_ss
tp = t(ip)

% tempo de subida 10-90% (interpolando entre amostras)
i10 = find(y >= 0.1*y_ss, 1);
i90 = find(y >= 0.9*y_ss, 1);
t10 = interp1(y(i10-1:i10), t(i10-1:i10), 0.1*y_ss);
t90 = interp1(y(i90-1:i90), t(i90-1:i90), 0.9*y_ss);
tr = t90 - t10

% tempo de acomodação (2%)
fora = find(abs(y - y_ss) > 0.02*y_ss);
ts = t(fora(end) + 1)

% esforço de controle (não vai a zero: a planta não tem integrador)
[umin, imin] = min(u);
u0   = u(1)
umin
u_ss = u(end)
var_total = max(u) - min(u)
du_max    = max(abs(diff(u)))

%% Gráficos
figure(1);
plot_y = stairs(t, y, 'b', 'LineWidth', 1.5); hold on;
plot_r = plot(t, r, 'k--');
faixa  = plot(t, 1.02*r, 'r:', t, 0.98*r, 'r:');
plot(tp, yp, 'ro', 'MarkerFaceColor', 'r');
text(tp + 0.1, yp, sprintf('M_p = %.1f%%, t_p = %.2f s', Mp, tp));
grid on; xlabel('t [s]'); ylabel('y[k]');
title('Resposta ao degrau em malha fechada');
legend([plot_y plot_r faixa(1)], {'y[k]', 'r[k]', 'faixa \pm2%'}, 'Location', 'southeast');

figure(2);
stairs(t, u, 'm', 'LineWidth', 1.5); hold on;
plot(t, 0*t, 'k--');
plot([t(1) t(imin)], [u(1) umin], 'ko', 'MarkerFaceColor', 'k');
text(t(1) + 0.1, u(1), sprintf('u(0) = %.2f', u(1)));
text(t(imin) + 0.1, umin, sprintf('u_{min} = %.2f em t = %.1f s', umin, t(imin)));
grid on; xlabel('t [s]'); ylabel('u[k]');
title('Esforço de controle');

figure(3);
pzmap(sysd, 'b', syscl, 'r'); zgrid; axis equal;
title('Polos e zeros: malha aberta (azul) e malha fechada (vermelho)');

%% Salvar figuras
pasta = fullfile('figuras', 'q2');
if ~exist(pasta, 'dir'), mkdir(pasta); end
saveas(figure(1), fullfile(pasta, 'resposta_degrau.png'));
saveas(figure(2), fullfile(pasta, 'sinal_controle.png'));
saveas(figure(3), fullfile(pasta, 'polos_z.png'));
