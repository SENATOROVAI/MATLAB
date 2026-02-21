clc;
clear;
close all;



s = tf('s');
Tu = 1;

% Передаточная функция разомкнутой системы
W_open = 1/(Tu*s*(Tu*s+1));

% Замыкаем систему по единичной обратной связи
W_closed = feedback(W_open,1);

% Строим переходной процесс
figure;
[y,t] = step(W_closed);
plot(t,y,'LineWidth',1.5);
grid on;
title('Переходной процесс — Линейный оптимум');
xlabel('Время, c');
ylabel('y(t)');


folder = fullfile(pwd,'figures');

if ~exist(folder,'dir')
    mkdir(folder);
end

saveas(gcf, fullfile(folder,'Fig_1_1_StepResponse_Linear.png'));

info2 = stepinfo(W_closed,'SettlingTimeThreshold',0.02);
info5 = stepinfo(W_closed,'SettlingTimeThreshold',0.05);

tp2 = info2.SettlingTime;
tp5 = info5.SettlingTime;
overshoot = info2.Overshoot;

fprintf('tp2 = %.4f c\n', tp2);
fprintf('tp5 = %.4f c\n', tp5);
fprintf('Перерегулирование = %.2f %%\n', overshoot);

% Биномиальный оптимум

W_bin = 1/(Tu*s + 1)^2;

figure;
[y_bin,t_bin] = step(W_bin);
plot(t_bin,y_bin,'LineWidth',1.5);
grid on;
title('Переходной процесс — Биномиальный оптимум');
xlabel('Время, c');
ylabel('y(t)');

saveas(gcf, fullfile(folder,'Fig_1_3_StepResponse_Binomial.png'));

info2_bin = stepinfo(W_bin,'SettlingTimeThreshold',0.02);
info5_bin = stepinfo(W_bin,'SettlingTimeThreshold',0.05);

tp2_bin = info2_bin.SettlingTime;
tp5_bin = info5_bin.SettlingTime;
overshoot_bin = info2_bin.Overshoot;

fprintf('\nБиномиальный оптимум\n');
fprintf('tp2 = %.4f c\n', tp2_bin);
fprintf('tp5 = %.4f c\n', tp5_bin);
fprintf('Перерегулирование = %.2f %%\n', overshoot_bin);

W_mod = 1/(2*Tu*s*(Tu*s+1));

W_mod_closed = feedback(W_mod,1);

figure;
[y_mod,t_mod] = step(W_mod_closed);
plot(t_mod,y_mod,'LineWidth',1.5);
grid on;
title('Переходной процесс — Оптимум по модулю');
xlabel('Время, c');
ylabel('y(t)');

saveas(gcf, fullfile(folder,'Fig_1_5_StepResponse_Modulus.png'));

info2_mod = stepinfo(W_mod_closed,'SettlingTimeThreshold',0.02);
info5_mod = stepinfo(W_mod_closed,'SettlingTimeThreshold',0.05);

tp2_mod = info2_mod.SettlingTime;
tp5_mod = info5_mod.SettlingTime;
overshoot_mod = info2_mod.Overshoot;

fprintf('\nОптимум по модулю\n');
fprintf('tp2 = %.4f c\n', tp2_mod);
fprintf('tp5 = %.4f c\n', tp5_mod);
fprintf('Перерегулирование = %.2f %%\n', overshoot_mod);

W_sym = 1/(Tu*s*(Tu*s+1)^2);
W_sym_closed = feedback(W_sym,1);

figure;
[y_sym,t_sym] = step(W_sym_closed);
plot(t_sym,y_sym,'LineWidth',1.5);
grid on;
title('Переходной процесс — Симметричный оптимум');
xlabel('Время, c');
ylabel('y(t)');

saveas(gcf, fullfile(folder,'Fig_1_6_StepResponse_Symmetric.png'));

info2_sym = stepinfo(W_sym_closed,'SettlingTimeThreshold',0.02);
info5_sym = stepinfo(W_sym_closed,'SettlingTimeThreshold',0.05);

tp2_sym = info2_sym.SettlingTime;
tp5_sym = info5_sym.SettlingTime;
overshoot_sym = info2_sym.Overshoot;

fprintf('\nСимметричный оптимум\n');
fprintf('tp2 = %.4f c\n', tp2_sym);
fprintf('tp5 = %.4f c\n', tp5_sym);
fprintf('Перерегулирование = %.2f %%\n', overshoot_sym);

W_ast = 1/(Tu*s)^3;
W_ast_closed = feedback(W_ast,1);

figure;
[y_ast,t_ast] = step(W_ast_closed);
plot(t_ast,y_ast,'LineWidth',1.5);
grid on;
title('Переходной процесс — Астатизм 3-го порядка');
xlabel('Время, c');
ylabel('y(t)');

saveas(gcf, fullfile(folder,'Fig_1_7_StepResponse_Astatism3.png'));

info2_ast = stepinfo(W_ast_closed,'SettlingTimeThreshold',0.02);
info5_ast = stepinfo(W_ast_closed,'SettlingTimeThreshold',0.05);

tp2_ast = info2_ast.SettlingTime;
tp5_ast = info5_ast.SettlingTime;
overshoot_ast = info2_ast.Overshoot;

fprintf('\nАстатизм 3 порядка\n');
fprintf('tp2 = %.4f c\n', tp2_ast);
fprintf('tp5 = %.4f c\n', tp5_ast);
fprintf('Перерегулирование = %.2f %%\n', overshoot_ast);

figure;
margin(W_open);
grid on;
title('ЛАЧХ и ФЧХ — Линейный оптимум');

saveas(gcf, fullfile(folder,'Fig_1_8_Bode_Linear.png'));

[gm, pm, wg, wp] = margin(W_open);

fprintf('\nЧастотные характеристики (Линейный оптимум)\n');
fprintf('Запас по амплитуде = %.2f dB\n', 20*log10(gm));
fprintf('Запас по фазе = %.2f град\n', pm);

t = 0:0.01:20;

v = 1;
a = 1;

g1 = v*t;        % g = vt
g2 = a*t.^2;     % g = at^2

% Линейный оптимум
y1_lin = lsim(W_closed,g1,t);
y2_lin = lsim(W_closed,g2,t);

figure;
plot(t,g1,'--',t,y1_lin,'LineWidth',1.5);
grid on;
title('Линейный оптимум — g = vt');
saveas(gcf, fullfile(folder,'Fig_Table2_Linear_vt.png'));

figure;
plot(t,g2,'--',t,y2_lin,'LineWidth',1.5);
grid on;
title('Линейный оптимум — g = at^2');
saveas(gcf, fullfile(folder,'Fig_Table2_Linear_at2.png'));




