%%  ��д�ߣ�
%%  ʱ�䣺 �� �� ��
%%  ����IEEE39ϵͳ�ĵ���ϵͳ��·����
clear;
clc;                                                                       %������ݱ����������д���
mpc = IEEE39;                                                              %��ȡ�����ļ� 
tic;                                                                       %��ʱ��ʼ       
f = input('���뽫������·���ϵĽڵ�f��');   
disp('1�������·��2�ǵ���ӵض�·��3�������·��4�������·�ӵ�');
type = input('��ѡ��������ͣ�');
branch = sparse(mpc.branch);                                               %��ȡĸ�����ݣ�������ϡ�����
generator = sparse(mpc.gen);                                               %��ȡ��������ݣ�������ϡ�����
a = -0.5+sqrt(3)*1i/2;                                                     %��ת���ӣ�a=e^(j120)
[b,c] = size(branch);                                                      %��ȡĸ�߾����������
nodes = 39;                                                                %IEEE39ϵͳ�Ľڵ���Ŀ39
%% �γ������ɾ���
bBij1 = -1./branch(:,5);                                                   %ĸ�ߵ�������
gBij1 = -1./generator(:,2);                                                %�������������
B1 = sparse(branch(:,3),branch(:,4),-bBij1,nodes,nodes);                   %�����ǷǶ���Ԫ�ػ�����
B1 = B1 + sparse(branch(:,4),branch(:,3),-bBij1,nodes,nodes);              %�����ǷǶ���Ԫ�ػ�����
B1 = B1 + sparse(branch(:,3),branch(:,3),bBij1,nodes,nodes);               %���Խ�Ԫ���Ե���
B1 = B1 + sparse(branch(:,4),branch(:,4),bBij1,nodes,nodes);
B1 = B1 + sparse(generator(:,1),generator(:,1),gBij1,nodes,nodes);
Y1 = 1i*B1;
%% �γɸ����ɾ���
bBij2 = -1./branch(:,5);                                                    %ĸ�ߵĸ�����
gBij2 = -1./generator(:,3);                                                 %������ĸ�����
B2 = sparse(branch(:,3),branch(:,4),-bBij2,nodes,nodes);
B2 = B2 + sparse(branch(:,4),branch(:,3),-bBij2,nodes,nodes);
B2 = B2 + sparse(branch(:,3),branch(:,3),bBij2,nodes,nodes);
B2 = B2 + sparse(branch(:,4),branch(:,4),bBij2,nodes,nodes);
B2 = B2 + sparse(generator(:,1),generator(:,1),gBij2,nodes,nodes);
Y2 = 1i*B2;
%% �γ������ɾ���
bBij0 = -1./branch(:,6);                                                    %ĸ�ߵ�������
B0 = sparse(branch(:,3),branch(:,4),-bBij0,nodes,nodes);
B0 = B0 + sparse(branch(:,4),branch(:,3),-bBij0,nodes,nodes);
B0 = B0 + sparse(branch(:,3),branch(:,3),bBij0,nodes,nodes);
B0 = B0 + sparse(branch(:,4),branch(:,4),bBij0,nodes,nodes);
%Yn,d��ѹ�������迹�����������ֵ���磬��������������ɾ��������
for i = 1:b                                                                 
    if (branch(i,2)== 3)
        B0(branch(i,3),branch(i,4)) = 0;                                   %Yn,d��ѹ��������ĸ�����ӻ�������0
        B0(branch(i,4),branch(i,3)) = 0;
        B0(branch(i,4),branch(i,4)) = 0;                                   %Yn,d��ѹ����ڵ��Ե�����0 
    end
end
%���Խ���Ԫ�ز���Ϊ0������Ϊ1
for j=1:nodes
    if(B0(j,j)==0)
        B0(j,j)=1;
    end
end
Y0 = 1i*B0;
%% ���ø���ڵ㵼�ɾ�������·��f�йصĽڵ��迹����ĵ�f��Ԫ��
If1 = zeros(nodes,1);                                                      %39��1��0��
If1(f) =1 ;                                                                %��·�ڵ���1
%����������������ȡU����������Z����+
Z1 = Y1\If1;
Z2 = Y2\If1;
Z0 = Y0\If1;
%% �����·
if(type==1)
    disp('---------------�����·---------------');
    %������ϵ������·����
    If_1 = 1/Z1(f);
    disp('���ϴ������·����:');
    disp(If_1);
    %�ǹ��ϵ��·��ѹ
    U_1=1-Z1*If_1;
    disp('�ǹ��Ͻڴ��ĵ�ѹΪ:');
    [d1,f1] = size(U_1);
    for i = 1:d1
        fprintf('  %d    %f\n',i,U_1(i));
    end
    %disp(U_1);
    %�����֧·�Ķ�·����
    U_ij1 = sparse (diag(U_1)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_1));
    I_ij1 = U_ij1.*Y1;
    disp('�����·�ǹ��ϴ���֧·�ĵ���Ϊ:')
    disp(sparse(I_ij1));    
end
%% ����ӵض�·
if(type==2)
    disp('---------------����ӵض�·---------------');
    %����ӵصĸ������ֵ���,���������СΪ
    If_2_1 = 1/(Z1(f)+Z2(f)+Z0(f));
    disp('����ӵض�·���ϴ�������������ֵ��ȣ���Ϊ:');                                                          
    disp(If_2_1);
    disp('���������ΪA�࣬����ӵض�·���ϴ���·����:');
    disp('A��');disp(3*If_2_1);
    %����ӵض�·���ϴ������ѹ
    Uf_2_1 = 1-Z1(f)*If_2_1;
    Uf_2_2 = 0-Z2(f)*If_2_1;
    Uf_2_0 = 0-Z0(f)*If_2_1;
    disp('���������򣬸��������ѹ:');
    disp(Uf_2_1);disp(Uf_2_2);disp(Uf_2_0);
    %���ϵ�A���ѹΪ�㣬�ɶԳƷ����������BC����ӵص�ѹ
    Uf_2_A = 0;
    Uf_2_B = a^2 * Uf_2_1 + a * Uf_2_2 + Uf_2_0;
    Uf_2_C = a * Uf_2_1 + a^2 * Uf_2_2 + Uf_2_0;
    Uf_2 = [Uf_2_A ; Uf_2_B ; Uf_2_C];                    
    disp('���ϴ�A��B��C��������ѹ�ֱ�Ϊ:' );
    disp(Uf_2);
   
    %����ǹ��ϴ����ڵ�ĵ�ѹ
    U_2_1 = 1-Z1*If_2_1;
    U_2_2 = 0-Z2*If_2_1;
    U_2_0 = 0-Z0*If_2_1;
    U_2_A = U_2_1 + U_2_2 + U_2_0;
    U_2_B = a^2 * U_2_1 + a * U_2_2 + U_2_0;
    U_2_C = a * U_2_1 + a^2 * U_2_2 + U_2_0;
    disp('����ӵض�·�ǹ��ϴ��ڵ�A���ѹ:');
    [d2A,f2A] = size(U_2_A);
    for i = 1:d2A
        fprintf('  %d    %f\n',i,U_2_A(i));
    end
    disp('����ӵض�·�ǹ��ϴ��ڵ�B���ѹ:');
    [d2B,f2B] = size(U_2_B);
    for i = 1:d2B
        fprintf('  %d    %f\n',i,U_2_B(i));
    end
    disp('����ӵض�·�ǹ��ϴ��ڵ�C���ѹ:');
    [d2C,f2C] = size(U_2_C);
    for i = 1:d2C
        fprintf('  %d    %f\n',i,U_2_C(i));
    end
    %�����֧·�Ķ�·����
    U_ij_21 = sparse (diag(U_2_1)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_2_1));
    U_ij_22 = sparse (diag(U_2_2)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_2_2));
    U_ij_20 = sparse (diag(U_2_0)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_2_0));
    I_ij_21 = U_ij_21.*Y1;
    I_ij_22 = U_ij_22.*Y2;
    I_ij_20 = U_ij_20.*Y0;
    I_ij_2_A = I_ij_21 + I_ij_22 + I_ij_20;
    I_ij_2_B = a^2 * I_ij_21 + a * I_ij_22 + I_ij_20;
    I_ij_2_C = a * I_ij_21 + a^2 * I_ij_22 + I_ij_20;
    disp('����ӵض�·�ǹ��ϴ���֧·A�����:');
    disp(I_ij_2_A);
    disp('����ӵض�·�ǹ��ϴ���֧·B�����:');
    disp(I_ij_2_B);
    disp('����ӵض�·�ǹ��ϴ���֧·C�����:');
    disp(I_ij_2_C);
end
%% �����·
if(type==3)
    disp('---------------�����·---------------');
    %���������·�������
    If_3_1 = 1/(Z1(f)+Z2(f));
    If_3_2 = -If_3_1;
    If_3_0 = 0;
    disp('�����·���������:');
    disp(If_3_1);
    disp('�����·�ĸ������:');
    disp(If_3_2);
    disp('�����·���������:');
    disp(If_3_0);
    disp('���������ΪB��C�࣬�������·������Ķ�·����Ϊ:');
    disp('B��');disp(-1i * sqrt(3) * If_3_1);
    disp('C��');disp(1i * sqrt(3) * If_3_1);
    %�����·�����ѹ
    Uf_3_1 = 1-Z1(f)*If_3_1;
    Uf_3_2 = 0-Z2(f)*If_3_2;
    Uf_3_0 = 0-Z0(f)*If_3_0;
    disp('���ϵ㴦A������򣬸��������ѹ:');
    disp(Uf_3_1);disp(Uf_3_2);disp(Uf_3_0);
    %��ΪA����B���·
    Uf_3_A = Uf_3_1 + Uf_3_2 ;
    Uf_3_B = a^2 * Uf_3_1 + a * Uf_3_2 ;
    Uf_3_C = a * Uf_3_1 + a^2 * Uf_3_2 ;
    Uf_3 = [Uf_3_A ; Uf_3_B ; Uf_3_C]; 
    disp('���ϴ�A��B��C�����ѹ�ֱ�Ϊ:' );
    disp(Uf_3);
    
    %����ǹ��ϴ����ڵ�ĵ�ѹ
    U_3_1 = 1-Z1*If_3_1;
    U_3_2 = 0-Z2*If_3_2;
    U_3_0 = 0-Z0*If_3_0;
    U_3_A = U_3_1 + U_3_2 + U_3_0;
    U_3_B = a^2 * U_3_1 + a * U_3_2 + U_3_0;
    U_3_C = a * U_3_1 + a^2 * U_3_2 + U_3_0;
    disp('�����·�ǹ��ϴ��ڵ�A���ѹ:');
    [d3A,f3A] = size(U_3_A);
    for i = 1:d3A
        fprintf('  %d    %f\n',i,U_3_A(i));
    end
    disp('�����·�ǹ��ϴ��ڵ�B���ѹ:');
    [d3B,f3B] = size(U_3_B);
    for i = 1:d3B
        fprintf('  %d    %f\n',i,U_3_B(i));
    end
    disp('�����·�ǹ��ϴ��ڵ�C���ѹ:');
    [d3C,f3C] = size(U_3_C);
    for i = 1:d3C
        fprintf('  %d    %f\n',i,U_3_C(i));
    end
    %�����֧·�Ķ�·����
    U_ij_31 = sparse (diag(U_3_1)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_3_1));
    U_ij_32 = sparse (diag(U_3_2)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_3_2));
    U_ij_30 = sparse (diag(U_3_0)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_3_0));
    I_ij_31 = U_ij_31.*Y1;
    I_ij_32 = U_ij_32.*Y2;
    I_ij_30 = U_ij_30.*Y0;
    I_ij_3_A = I_ij_31 + I_ij_32 + I_ij_30;
    I_ij_3_B = a^2 * I_ij_31 + a * I_ij_32 + I_ij_30;
    I_ij_3_C = a * I_ij_31 + a^2 * I_ij_32 + I_ij_30;
    disp('�����·�ǹ��ϴ���֧·A�����:');
    disp(I_ij_3_A);
    disp('�����·�ǹ��ϴ���֧·B�����:');
    disp(I_ij_3_B);
    disp('�����·�ǹ��ϴ���֧·C�����:');
    disp(I_ij_3_C);
end
%% �����·�ӵ�
if(type==4)
    disp('---------------�����·�ӵ�---------------');
    %���������·�ӵظ������
    If_4_1 = 1/(Z1(f)+(Z2(f)*Z0(f))/(Z2(f)+Z0(f)));
    If_4_2 = -If_4_1 * Z0(f)/(Z2(f)+Z0(f));
    If_4_0 = -If_4_1 * Z2(f)/(Z2(f)+Z0(f));
    disp('�����·�ӵص��������:');
    disp(If_4_1);
    disp('�����·�ӵصĸ������:');
    disp(If_4_2);
    disp('�����·�ӵص��������:');
    disp(If_4_0);
    %B��C����ӵض�·����
    If_4 =If_4_1 * sqrt(3) * sqrt(1-Z2(f)*Z0(f)/( Z2(f) + Z0(f))^2);
    disp('���������ΪB��C��,�������·�ӵع�����Ķ�·����:');
    disp('B��');disp(If_4);
    disp('C��');disp(If_4);
    Uf_4_1 = 1-Z1(f)*If_4_1;
    Uf_4_2 = 0-Z2(f)*If_4_2;
    Uf_4_0 = 0-Z0(f)*If_4_0;
    disp('���ϵ㴦A������򣬸��������ѹ:');
    disp(Uf_4_1);disp(Uf_4_2);disp(Uf_4_0);
    %����ʵ�ʵ�ѹ
    Uf_4_A = Uf_4_1 + Uf_4_2 + Uf_4_0;
    Uf_4_B = a^2 * Uf_4_1 + a * Uf_4_2 + Uf_4_0;
    Uf_4_C = a * Uf_4_1 + a^2 * Uf_4_2 + Uf_4_0;
    Uf_4 = [Uf_4_A ; Uf_4_B ; Uf_4_C];             
    disp('���ϴ�A��B��C��������ѹ�ֱ�Ϊ:' );
    disp(Uf_4);
    
    %����ǹ��ϴ����ڵ�ĵ�ѹ
    U_4_1 = 1-Z1*If_4_1;
    U_4_2 = 0-Z2*If_4_2;
    U_4_0 = 0-Z0*If_4_0;
    U_4_A = U_4_1 + U_4_2 + U_4_0;
    U_4_B = a^2 * U_4_1 + a * U_4_2 + U_4_0;
    U_4_C = a * U_4_1 + a^2 * U_4_2 + U_4_0;
    disp('�����·�ӵطǹ��ϴ��ڵ�A���ѹ:');
    [d4A,f4A] = size(U_4_A);
    for i = 1:d4A
        fprintf('  %d    %f\n',i,U_4_A(i));
    end
    disp('�����·�ӵطǹ��ϴ��ڵ�B���ѹ:');
    [d4B,f4B] = size(U_4_B);
    for i = 1:d4B
        fprintf('  %d    %f\n',i,U_4_B(i));
    end
    disp('�����·�ӵطǹ��ϴ��ڵ�C���ѹ:');
    [d4C,f4C] = size(U_4_C);
    for i = 1:d4C
        fprintf('  %d    %f\n',i,U_4_C(i));
    end
    %�����֧·�Ķ�·����
    U_ij_41 = sparse (diag(U_4_1)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_4_1));
    U_ij_42 = sparse (diag(U_4_2)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_4_2));
    U_ij_40 = sparse (diag(U_4_0)*ones(nodes,nodes)-ones(nodes,nodes)*diag(U_4_0));
    I_ij_41 = U_ij_41.*Y1;
    I_ij_42 = U_ij_42.*Y2;
    I_ij_40 = U_ij_40.*Y0;
    I_ij_4_A = I_ij_41 + I_ij_42 + I_ij_40;
    I_ij_4_B = a^2 * I_ij_41 + a * I_ij_42 + I_ij_40;
    I_ij_4_C = a * I_ij_41 + a^2 * I_ij_42 + I_ij_40;
    disp('�����·�ӵطǹ��ϴ���֧·A�����:');
    disp(I_ij_4_A);
    disp('�����·�ӵطǹ��ϴ���֧·B�����:');
    disp(I_ij_4_B);
    disp('�����·�ӵطǹ��ϴ���֧·C�����:');
    disp(I_ij_4_C);
end 
toc;