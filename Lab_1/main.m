clc;
clear;
close all;

%% Параметры
s = tf('s');
Tu = 1;

folder = fullfile(pwd,'figures');
if ~exist(folder,'dir')
    mkdir(folder);
end

%% ============================================
%% 1. ЛИНЕЙНЫЙ ОПТИМУМ
%% ============================================

W_open_lin = 1/(Tu*s*(Tu*s+1));
W_closed_lin = feedback(W_open_lin,1);

%% Переходный процесс
figure;
step(W_closed_lin);
grid on;
title('Линейный оптимум');
saveas(gcf, fullfile(folder,'Lin_Step.png'));

info2 = stepinfo(W_closed_lin,'SettlingTimeThreshold',0.02);
info5 = stepinfo(W_closed_lin,'SettlingTimeThreshold',0.05);

%% ============================================
%% 2. БИНОМИАЛЬНЫЙ ОПТИМУМ (ЭТАЛОННАЯ ЗАМКНУТАЯ МОДЕЛЬ)
%% ============================================

W_closed_bin = 1/(Tu*s + 1)^2;

% Восстанавливаем разомкнутую систему для построения ЛАЧХ
W_open_bin = W_closed_bin/(1 - W_closed_bin);

figure;
step(W_closed_bin);
grid on;
title('Биномиальный оптимум');
saveas(gcf, fullfile(folder,'Bin_Step.png'));

info2_bin = stepinfo(W_closed_bin,'SettlingTimeThreshold',0.02);
info5_bin = stepinfo(W_closed_bin,'SettlingTimeThreshold',0.05);

%% ============================================
%% 3. ОПТИМУМ ПО МОДУЛЮ
%% ============================================

W_open_mod = 1/(2*Tu*s*(Tu*s+1));
W_closed_mod = feedback(W_open_mod,1);

figure;
step(W_closed_mod);
grid on;
title('Оптимум по модулю');
saveas(gcf, fullfile(folder,'Mod_Step.png'));

info2_mod = stepinfo(W_closed_mod,'SettlingTimeThreshold',0.02);
info5_mod = stepinfo(W_closed_mod,'SettlingTimeThreshold',0.05);

%% ============================================
%% 4. СИММЕТРИЧНЫЙ ОПТИМУМ
%% ============================================

W_open_sym = 1/(Tu*s*(Tu*s+1)^2);
W_closed_sym = feedback(W_open_sym,1);

figure;
step(W_closed_sym);
grid on;
title('Симметричный оптимум');
saveas(gcf, fullfile(folder,'Sym_Step.png'));

info2_sym = stepinfo(W_closed_sym,'SettlingTimeThreshold',0.02);
info5_sym = stepinfo(W_closed_sym,'SettlingTimeThreshold',0.05);

%% ============================================
%% 5. АСТАТИЗМ 3 ПОРЯДКА
%% ============================================

W_open_ast = 1/(Tu*s)^3;
W_closed_ast = feedback(W_open_ast,1);

figure;
step(W_closed_ast);
grid on;
title('Астатизм 3 порядка');
saveas(gcf, fullfile(folder,'Ast_Step.png'));

info2_ast = stepinfo(W_closed_ast,'SettlingTimeThreshold',0.02);
info5_ast = stepinfo(W_closed_ast,'SettlingTimeThreshold',0.05);

%% ============================================
%% ТАБЛИЦА 1
%% ============================================

fprintf('\nТАБЛИЦА 1\n');
fprintf('Модель | tp5/Tu | tp2/Tu | Overshoot\n');

fprintf('Линейный | %.4f | %.4f | %.2f\n', ...
    info5.SettlingTime/Tu, info2.SettlingTime/Tu, info2.Overshoot);

fprintf('Биномиальный | %.4f | %.4f | %.2f\n', ...
    info5_bin.SettlingTime/Tu, info2_bin.SettlingTime/Tu, info2_bin.Overshoot);

fprintf('По модулю | %.4f | %.4f | %.2f\n', ...
    info5_mod.SettlingTime/Tu, info2_mod.SettlingTime/Tu, info2_mod.Overshoot);

fprintf('Симметричный | %.4f | %.4f | %.2f\n', ...
    info5_sym.SettlingTime/Tu, info2_sym.SettlingTime/Tu, info2_sym.Overshoot);

fprintf('Астатизм 3 | %.4f | %.4f | %.2f\n', ...
    info5_ast.SettlingTime/Tu, info2_ast.SettlingTime/Tu, info2_ast.Overshoot);

%% ============================================
%% ЛАЧХ И ЗАПАСЫ УСТОЙЧИВОСТИ
%% ============================================

models = {W_open_lin, W_open_bin, W_open_mod, W_open_sym, W_open_ast};
names = {'Lin','Bin','Mod','Sym','Ast'};

fprintf('\nТАБЛИЦА 3\n');
fprintf('Модель | GM(dB) | PM(deg) | M\n');

for i = 1:length(models)
    
    figure;
    margin(models{i});
    grid on;
    saveas(gcf, fullfile(folder, strcat(names{i},'_Bode.png')));
    
    [gm,pm,~,~] = margin(models{i});
    
    Wcl = feedback(models{i},1);
    [mag,~,~] = bode(Wcl);
    M = max(squeeze(mag));
    
    fprintf('%s | %.2f | %.2f | %.2f\n', ...
        names{i}, 20*log10(gm), pm, M);
end

%% ============================================
%% ТАБЛИЦА 2 (g = vt, g = at^2)
%% ============================================

t = 0:0.01:20;
v = 1;
a = 1;

g1 = v*t;
g2 = a*t.^2;

systems = {W_closed_lin, W_closed_bin, W_closed_mod, W_closed_sym, W_closed_ast};

fprintf('\nТАБЛИЦА 2 (ошибки установившиеся)\n');

for i = 1:length(systems)
    
    y1 = lsim(systems{i},g1,t);
    y2 = lsim(systems{i},g2,t);
    
    e1 = abs(g1(end)-y1(end));
    e2 = abs(g2(end)-y2(end));
    
    fprintf('%s | e(vt)=%.4f | e(at2)=%.4f\n', names{i}, e1, e2);
end