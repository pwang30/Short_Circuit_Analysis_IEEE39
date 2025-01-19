%%  编写者：
%%  时间： 年 月 日
%%  基于IEEE39系统的电力系统短路计算
clear;
clc;                                                                       %清空数据变量及命令行窗口
mpc = IEEE39;                                                              %读取数据文件 
tic;                                                                       %计时开始       
f = input('输入将发生短路故障的节点f：');   
disp('1是三相短路，2是单相接地短路，3是两相短路，4是两相短路接地');
type = input('请选择故障类型：');
branch = sparse(mpc.branch);                                               %获取母线数据，并建立稀疏矩阵
generator = sparse(mpc.gen);                                               %获取发电机数据，并建立稀疏矩阵
a = -0.5+sqrt(3)*1i/2;                                                     %旋转因子，a=e^(j120)
[b,c] = size(branch);                                                      %获取母线矩阵的行列数
nodes = 39;                                                                %IEEE39系统的节点数目39
%% 形成正序导纳矩阵
bBij1 = -1./branch(:,5);                                                   %母线的正序导纳
gBij1 = -1./generator(:,2);                                                %发电机的正序导纳
B1 = sparse(branch(:,3),branch(:,4),-bBij1,nodes,nodes);                   %上三角非对象元素互导纳
B1 = B1 + sparse(branch(:,4),branch(:,3),-bBij1,nodes,nodes);              %上三角非对象元素互导纳
B1 = B1 + sparse(branch(:,3),branch(:,3),bBij1,nodes,nodes);               %主对角元素自导纳
B1 = B1 + sparse(branch(:,4),branch(:,4),bBij1,nodes,nodes);
B1 = B1 + sparse(generator(:,1),generator(:,1),gBij1,nodes,nodes);
Y1 = 1i*B1;
%% 形成负序导纳矩阵
bBij2 = -1./branch(:,5);                                                    %母线的负序导纳
gBij2 = -1./generator(:,3);                                                 %发电机的负序导纳
B2 = sparse(branch(:,3),branch(:,4),-bBij2,nodes,nodes);
B2 = B2 + sparse(branch(:,4),branch(:,3),-bBij2,nodes,nodes);
B2 = B2 + sparse(branch(:,3),branch(:,3),bBij2,nodes,nodes);
B2 = B2 + sparse(branch(:,4),branch(:,4),bBij2,nodes,nodes);
B2 = B2 + sparse(generator(:,1),generator(:,1),gBij2,nodes,nodes);
Y2 = 1i*B2;
%% 形成零序导纳矩阵
bBij0 = -1./branch(:,6);                                                    %母线的零序导纳
B0 = sparse(branch(:,3),branch(:,4),-bBij0,nodes,nodes);
B0 = B0 + sparse(branch(:,4),branch(:,3),-bBij0,nodes,nodes);
B0 = B0 + sparse(branch(:,3),branch(:,3),bBij0,nodes,nodes);
B0 = B0 + sparse(branch(:,4),branch(:,4),bBij0,nodes,nodes);
%Yn,d变压器零序阻抗不计入零序等值网络，将其参数从零序导纳矩阵中清除
for i = 1:b                                                                 
    if (branch(i,2)== 3)
        B0(branch(i,3),branch(i,4)) = 0;                                   %Yn,d变压器与其他母线连接互导纳置0
        B0(branch(i,4),branch(i,3)) = 0;
        B0(branch(i,4),branch(i,4)) = 0;                                   %Yn,d变压器侧节点自导纳置0 
    end
end
%主对角线元素不能为0，修正为1
for j=1:nodes
    if(B0(j,j)==0)
        B0(j,j)=1;
    end
end
Y0 = 1i*B0;
%% 利用各序节点导纳矩阵计算短路点f有关的节点阻抗矩阵的第f列元素
If1 = zeros(nodes,1);                                                      %39行1列0阵
If1(f) =1 ;                                                                %短路节点置1
%运用左除处理矩阵，求取U向量，即得Z向量+
Z1 = Y1\If1;
Z2 = Y2\If1;
Z0 = Y0\If1;
%% 三相短路
if(type==1)
    disp('---------------三相短路---------------');
    %计算故障点三相短路电流
    If_1 = 1/Z1(f);
    disp('故障处三相短路电流:');
    disp(If_1);
    %非故障点短路电压
    U_1=1-Z1*If_1;
    disp('非故障节处的电压为:');
    [d1,f1] = size(U_1);
    for i = 1:d1
        fprintf('  %d    %f\n',i,U_1(i));
    end
    %disp(U_1);
    %计算各支路的短路电流
    U_ij1 = sparse (diag(U_1)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_1));
    I_ij1 = U_ij1.*Y1;
    disp('三相短路非故障处各支路的电流为:')
    disp(sparse(I_ij1));    
end
%% 单相接地短路
if(type==2)
    disp('---------------单相接地短路---------------');
    %单相接地的各序电流值相等,正序电流大小为
    If_2_1 = 1/(Z1(f)+Z2(f)+Z0(f));
    disp('单相接地短路故障处正负零各序电流值相等，均为:');                                                          
    disp(If_2_1);
    disp('假设故障相为A相，单相接地短路故障处短路电流:');
    disp('A相');disp(3*If_2_1);
    %单相接地短路故障处各序电压
    Uf_2_1 = 1-Z1(f)*If_2_1;
    Uf_2_2 = 0-Z2(f)*If_2_1;
    Uf_2_0 = 0-Z0(f)*If_2_1;
    disp('故障相正序，负序，零序电压:');
    disp(Uf_2_1);disp(Uf_2_2);disp(Uf_2_0);
    %故障点A相电压为零，由对称分量法并求出BC单相接地电压
    Uf_2_A = 0;
    Uf_2_B = a^2 * Uf_2_1 + a * Uf_2_2 + Uf_2_0;
    Uf_2_C = a * Uf_2_1 + a^2 * Uf_2_2 + Uf_2_0;
    Uf_2 = [Uf_2_A ; Uf_2_B ; Uf_2_C];                    
    disp('故障处A、B、C三相各相电压分别为:' );
    disp(Uf_2);
   
    %计算非故障处各节点的电压
    U_2_1 = 1-Z1*If_2_1;
    U_2_2 = 0-Z2*If_2_1;
    U_2_0 = 0-Z0*If_2_1;
    U_2_A = U_2_1 + U_2_2 + U_2_0;
    U_2_B = a^2 * U_2_1 + a * U_2_2 + U_2_0;
    U_2_C = a * U_2_1 + a^2 * U_2_2 + U_2_0;
    disp('单相接地短路非故障处节点A相电压:');
    [d2A,f2A] = size(U_2_A);
    for i = 1:d2A
        fprintf('  %d    %f\n',i,U_2_A(i));
    end
    disp('单相接地短路非故障处节点B相电压:');
    [d2B,f2B] = size(U_2_B);
    for i = 1:d2B
        fprintf('  %d    %f\n',i,U_2_B(i));
    end
    disp('单相接地短路非故障处节点C相电压:');
    [d2C,f2C] = size(U_2_C);
    for i = 1:d2C
        fprintf('  %d    %f\n',i,U_2_C(i));
    end
    %计算各支路的短路电流
    U_ij_21 = sparse (diag(U_2_1)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_2_1));
    U_ij_22 = sparse (diag(U_2_2)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_2_2));
    U_ij_20 = sparse (diag(U_2_0)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_2_0));
    I_ij_21 = U_ij_21.*Y1;
    I_ij_22 = U_ij_22.*Y2;
    I_ij_20 = U_ij_20.*Y0;
    I_ij_2_A = I_ij_21 + I_ij_22 + I_ij_20;
    I_ij_2_B = a^2 * I_ij_21 + a * I_ij_22 + I_ij_20;
    I_ij_2_C = a * I_ij_21 + a^2 * I_ij_22 + I_ij_20;
    disp('单相接地短路非故障处各支路A相电流:');
    disp(I_ij_2_A);
    disp('单相接地短路非故障处各支路B相电流:');
    disp(I_ij_2_B);
    disp('单相接地短路非故障处各支路C相电流:');
    disp(I_ij_2_C);
end
%% 两相短路
if(type==3)
    disp('---------------两相短路---------------');
    %计算两相短路各序电流
    If_3_1 = 1/(Z1(f)+Z2(f));
    If_3_2 = -If_3_1;
    If_3_0 = 0;
    disp('两相短路的正序电流:');
    disp(If_3_1);
    disp('两相短路的负序电流:');
    disp(If_3_2);
    disp('两相短路的零序电流:');
    disp(If_3_0);
    disp('假设故障相为B、C相，则两相短路故障相的短路电流为:');
    disp('B相');disp(-1i * sqrt(3) * If_3_1);
    disp('C相');disp(1i * sqrt(3) * If_3_1);
    %两相短路各序电压
    Uf_3_1 = 1-Z1(f)*If_3_1;
    Uf_3_2 = 0-Z2(f)*If_3_2;
    Uf_3_0 = 0-Z0(f)*If_3_0;
    disp('故障点处A相的正序，负序，零序电压:');
    disp(Uf_3_1);disp(Uf_3_2);disp(Uf_3_0);
    %设为A相与B相短路
    Uf_3_A = Uf_3_1 + Uf_3_2 ;
    Uf_3_B = a^2 * Uf_3_1 + a * Uf_3_2 ;
    Uf_3_C = a * Uf_3_1 + a^2 * Uf_3_2 ;
    Uf_3 = [Uf_3_A ; Uf_3_B ; Uf_3_C]; 
    disp('故障处A、B、C各相电压分别为:' );
    disp(Uf_3);
    
    %计算非故障处各节点的电压
    U_3_1 = 1-Z1*If_3_1;
    U_3_2 = 0-Z2*If_3_2;
    U_3_0 = 0-Z0*If_3_0;
    U_3_A = U_3_1 + U_3_2 + U_3_0;
    U_3_B = a^2 * U_3_1 + a * U_3_2 + U_3_0;
    U_3_C = a * U_3_1 + a^2 * U_3_2 + U_3_0;
    disp('两相短路非故障处节点A相电压:');
    [d3A,f3A] = size(U_3_A);
    for i = 1:d3A
        fprintf('  %d    %f\n',i,U_3_A(i));
    end
    disp('两相短路非故障处节点B相电压:');
    [d3B,f3B] = size(U_3_B);
    for i = 1:d3B
        fprintf('  %d    %f\n',i,U_3_B(i));
    end
    disp('两相短路非故障处节点C相电压:');
    [d3C,f3C] = size(U_3_C);
    for i = 1:d3C
        fprintf('  %d    %f\n',i,U_3_C(i));
    end
    %计算各支路的短路电流
    U_ij_31 = sparse (diag(U_3_1)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_3_1));
    U_ij_32 = sparse (diag(U_3_2)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_3_2));
    U_ij_30 = sparse (diag(U_3_0)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_3_0));
    I_ij_31 = U_ij_31.*Y1;
    I_ij_32 = U_ij_32.*Y2;
    I_ij_30 = U_ij_30.*Y0;
    I_ij_3_A = I_ij_31 + I_ij_32 + I_ij_30;
    I_ij_3_B = a^2 * I_ij_31 + a * I_ij_32 + I_ij_30;
    I_ij_3_C = a * I_ij_31 + a^2 * I_ij_32 + I_ij_30;
    disp('两相短路非故障处各支路A相电流:');
    disp(I_ij_3_A);
    disp('两相短路非故障处各支路B相电流:');
    disp(I_ij_3_B);
    disp('两相短路非故障处各支路C相电流:');
    disp(I_ij_3_C);
end
%% 两相短路接地
if(type==4)
    disp('---------------两相短路接地---------------');
    %计算两相短路接地各序电流
    If_4_1 = 1/(Z1(f)+(Z2(f)*Z0(f))/(Z2(f)+Z0(f)));
    If_4_2 = -If_4_1 * Z0(f)/(Z2(f)+Z0(f));
    If_4_0 = -If_4_1 * Z2(f)/(Z2(f)+Z0(f));
    disp('两相短路接地的正序电流:');
    disp(If_4_1);
    disp('两相短路接地的负序电流:');
    disp(If_4_2);
    disp('两相短路接地的零序电流:');
    disp(If_4_0);
    %B、C两相接地短路电流
    If_4 =If_4_1 * sqrt(3) * sqrt(1-Z2(f)*Z0(f)/( Z2(f) + Z0(f))^2);
    disp('假设故障相为B、C相,则两相短路接地故障相的短路电流:');
    disp('B相');disp(If_4);
    disp('C相');disp(If_4);
    Uf_4_1 = 1-Z1(f)*If_4_1;
    Uf_4_2 = 0-Z2(f)*If_4_2;
    Uf_4_0 = 0-Z0(f)*If_4_0;
    disp('故障点处A相的正序，负序，零序电压:');
    disp(Uf_4_1);disp(Uf_4_2);disp(Uf_4_0);
    %各相实际电压
    Uf_4_A = Uf_4_1 + Uf_4_2 + Uf_4_0;
    Uf_4_B = a^2 * Uf_4_1 + a * Uf_4_2 + Uf_4_0;
    Uf_4_C = a * Uf_4_1 + a^2 * Uf_4_2 + Uf_4_0;
    Uf_4 = [Uf_4_A ; Uf_4_B ; Uf_4_C];             
    disp('故障处A、B、C三相各相电压分别为:' );
    disp(Uf_4);
    
    %计算非故障处各节点的电压
    U_4_1 = 1-Z1*If_4_1;
    U_4_2 = 0-Z2*If_4_2;
    U_4_0 = 0-Z0*If_4_0;
    U_4_A = U_4_1 + U_4_2 + U_4_0;
    U_4_B = a^2 * U_4_1 + a * U_4_2 + U_4_0;
    U_4_C = a * U_4_1 + a^2 * U_4_2 + U_4_0;
    disp('两相短路接地非故障处节点A相电压:');
    [d4A,f4A] = size(U_4_A);
    for i = 1:d4A
        fprintf('  %d    %f\n',i,U_4_A(i));
    end
    disp('两相短路接地非故障处节点B相电压:');
    [d4B,f4B] = size(U_4_B);
    for i = 1:d4B
        fprintf('  %d    %f\n',i,U_4_B(i));
    end
    disp('两相短路接地非故障处节点C相电压:');
    [d4C,f4C] = size(U_4_C);
    for i = 1:d4C
        fprintf('  %d    %f\n',i,U_4_C(i));
    end
    %计算各支路的短路电流
    U_ij_41 = sparse (diag(U_4_1)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_4_1));
    U_ij_42 = sparse (diag(U_4_2)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_4_2));
    U_ij_40 = sparse (diag(U_4_0)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_4_0));
    I_ij_41 = U_ij_41.*Y1;
    I_ij_42 = U_ij_42.*Y2;
    I_ij_40 = U_ij_40.*Y0;
    I_ij_4_A = I_ij_41 + I_ij_42 + I_ij_40;
    I_ij_4_B = a^2 * I_ij_41 + a * I_ij_42 + I_ij_40;
    I_ij_4_C = a * I_ij_41 + a^2 * I_ij_42 + I_ij_40;
    disp('两相短路接地非故障处各支路A相电流:');
    disp(I_ij_4_A);
    disp('两相短路接地非故障处各支路B相电流:');
    disp(I_ij_4_B);
    disp('两相短路接地非故障处各支路C相电流:');
    disp(I_ij_4_C);
end 
toc;