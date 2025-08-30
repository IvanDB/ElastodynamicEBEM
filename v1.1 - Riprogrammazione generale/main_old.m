close all
clear 
clc

format long e

warning off

%Definiamo il percorso delle cartelle in cui Matlab dovrà ricercare
%i file da leggere e le funzioni da richiamare

addpath ./old_function
addpath ./old_function/src ./old_function/src/src_input ./old_function/src/src_quadrature ./old_function/src/src_messages
addpath ./old_function/src/src_rhs
addpath ./old_function/Post_processing
addpath ./old_function/input_functions
addpath ./old_function/quadrature

%% LETTURA dei PARAMETRI di INPUT

%Lettura del file di input
pb_param = read_input_file;

%La function read_input_file() restituisce in OUTPUT nella struct 
%pb_param i parametri del problema che sono memorizzati all'interno del
%file 'input_file'

%Verifica che il metodo risolutivo sia stato implementato
[err_flag,message] = check_implementation(pb_param);
if err_flag == 0
    error(message)
    return
end
%La function check_implementation() prende in INPUT la struct pb.param e 
%controlla che il metodo risolutivo scelto dall'utente sia stato 
%implementato. In particolare essa restitusice in OUTPUT:
%
%-la variabile err_flag che vale 0 se sono presenti alcuni errori.
% In caso contrario essa assume valore 1
%
%-la variabile message, che contiene un messaggio di errore nel caso in cui
% sia presente un errore. In caso contrario essa è una stringa vuota.


% %% CREAZIONE del FILE di OUTPUT
% 
% %Nome del file di output
% file_out=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_',...
%     num2str(kappa),'_',pb_param.module_type,...
%     '_',num2str(pb_param.deg_k),'_out','.txt');
% %Indice identificativo del file di output
% fid = fopen (file_out, 'w');
% 
% %Stampa a video i parametri del problema
% init_message(fid,pb_param.domain_phys,pb_param.domain_type,pb_param.lev,...
%     pb_param.module_type,pb_param.deg_k,kappa);


%% LETTURA della MESH SPAZIALE

%Lettura del file contenente la MESH SPAZIALE
domainMesh = time3D_read_space_mesh(pb_param.domain_type,pb_param.lev);

%La funzione time3D_read_space_mesh() prende in INPUT:
% -la variabile pb_param.domain_type contenente il TIPO DI DOMINIO
% -la variabile pb_param.lev contenente il LIVELLO DI RAFFINAMENTO
%
%e restituisce in OUTPUT la struct domainMesh


%% NODI e PESI per FORMULE di QUADRATURA 

%Costruiamo i NODI e i PESI della FORMULA di QUADRATURA

mxghp = 28;   %Massimo numero di nodi e pesi di quadratura

%Matrici contenenti i nodi e i pesi per formule di quadratura di Gauss-Hammer
[gha,ghw] = gausshammer(mxghp);

%% TIME-MARCHING

%Plot della discretizzazione del dominio considerato
%time3D_plot(domainMesh)

%Calcoliamo la funzione densità
density = time3D_history(pb_param,domainMesh,gha,ghw);


%% Plot della FUNZIONE DENSITA'
X = [(domainMesh.coordinates(domainMesh.triangles(:,1),1))';...
     (domainMesh.coordinates(domainMesh.triangles(:,2),1))';...
     (domainMesh.coordinates(domainMesh.triangles(:,3),1))'];
 
Y = [(domainMesh.coordinates(domainMesh.triangles(:,1),2))';...
     (domainMesh.coordinates(domainMesh.triangles(:,2),2))';...
     (domainMesh.coordinates(domainMesh.triangles(:,3),2))']; 
 
Z = [(domainMesh.coordinates(domainMesh.triangles(:,1),3))';...
     (domainMesh.coordinates(domainMesh.triangles(:,2),3))';...
     (domainMesh.coordinates(domainMesh.triangles(:,3),3))']; 
 
%Plot funzione densità sulla mesh relativamente all'istante 
%temporale di indice hk
hk=pb_param.Nt;
for j=1:3
    figure(j)
    fill3(X,Y,Z,density(j:3:end,hk))
    xlabel('Asse $x_1$','interpreter','latex','fontsize',14);
    ylabel('Asse $x_2$','interpreter','latex','fontsize',14);
    zlabel('Asse $x_3$','interpreter','latex','fontsize',14);
    title(strcat('$\phi_',num2str(j),'(x_1,x_2,x_3,T)$ con $ T=',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',14); 
    colormap(jet)
    colorbar
    daspect([1 1 1]);
    %caxis([0 3])
    %view(2)
    name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_phi',num2str(j),'_intervallo_temporale_',num2str(hk));
    print(name,'-dtiff')
    savefig(strcat(name,'.fig'))
    print(name,'-depsc2','-r500')
end

phi1=density(1:3:end,end);
phi1_min=min(phi1);
phi1_max=max(phi1);
max_diff1=max(abs(phi1-ones(length(phi1),1)));

phi2=density(2:3:end,end);
phi2_min=min(phi2);
phi2_max=max(phi2);
max_diff2=max(abs(phi2-ones(length(phi2),1)));

phi3=density(3:3:end,end);
phi3_min=min(phi3);
phi3_max=max(phi3);
max_diff3=max(abs(phi3-ones(length(phi3),1)));

name=strcat('Workspace_',pb_param.domain_type,'_',num2str(pb_param.lev),'_Nodi12_T',num2str(pb_param.T_fin),'_N',num2str(pb_param.Nt),'_Delta1_Beta1');
save(name)

% %% Densità sulla barretta 
% T=pb_param.T_fin;
% IND=[96 80 48 27 59 113];
% %IND=[608 625 571 539 560 592];
% for j=1:length(IND)
% ind=IND(j);
% for i=1:3
%     alpha = density(3*ind-(3-i),:);
%     tt = linspace(0,pb_param.T_fin,pb_param.Nt+1);
%     XX = linspace(0,pb_param.T_fin,2001);
%     YY=zeros(1,length(XX));
%     for l=1:pb_param.Nt-1
%         chi_l=@(t)((t >= tt(l) & t<tt(l+1)).*1);
%         YY=YY+chi_l(XX)*alpha(l);
%     end
%     chi_l=@(t)((t >= tt(pb_param.Nt) & t<=tt(pb_param.Nt+1)).*1);
%     YY=YY+chi_l(XX)*alpha(pb_param.Nt);
%     figure(i)
%     box off
%     plot(XX,YY)
%     xlim([0 T])
%     xlabel('Asse $t$','interpreter','latex','fontsize',13);
%     ylabel(strcat('$\phi_',num2str(i)','(x_1,x_2,x_3,t)$'),'interpreter','latex','fontsize',16);
%     title(strcat('Grafico di $\phi_',num2str(i),'(x_1,x_2,x_3,t)$ per $t\in[0,',num2str(T),']$ sul triangolo $',num2str(ind),'$'),'interpreter','latex','fontsize',12); 
%     name=strcat('phi',num2str(i),'_punto_triangolo_',num2str(ind));
%     print(name,'-dtiff')
%     savefig(strcat(name,'.fig'))
%     %close all
% end
% end
% 
% 
% ind=592;
% i=3;
% T=pb_param.T_fin;
% alpha = density(3*ind-(3-i),:);
%     tt = linspace(0,pb_param.T_fin,pb_param.Nt+1);
%     XX = linspace(0,pb_param.T_fin,2001);
%     YY=zeros(1,length(XX));
%     for l=1:pb_param.Nt-1
%         chi_l=@(t)((t >= tt(l) & t<tt(l+1)).*1);
%         YY=YY+chi_l(XX)*alpha(l);
%     end
%     chi_l=@(t)((t >= tt(pb_param.Nt) & t<=tt(pb_param.Nt+1)).*1);
%     YY=YY+chi_l(XX)*alpha(pb_param.Nt);
%     figure(i)
%     box off
%     plot(XX,YY,'linewidth',0.8)
%     xlim([0 T])
%     box off
%     set(gca,'TickLabelInterpreter','latex')
% %     xlabel('Asse $t$','interpreter','latex','fontsize',13);
% %     ylabel(strcat('$\phi_',num2str(i)','(x_1,x_2,x_3,t)$'),'interpreter','latex','fontsize',16);
% %     title(strcat('Grafico di $\phi_',num2str(i),'(x_1,x_2,x_3,t)$ per $t\in[0,',num2str(T),']$ sul triangolo $',num2str(ind),'$'),'interpreter','latex','fontsize',12); 
%     name=strcat('phi',num2str(i),'_punto_triangolo_',num2str(ind));
%     %print(name,'-dtiff')
%     savefig(strcat(name,'.fig'))
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% screenheiko1
%------------------------------------------------------------------
%------------------------------------------------------------------
%Comandi per plottare soluzione relativa a screenheiko1_0
% ind = [7 3 4];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% %xx = linspace(-0.5,0.5,length(ind));
% xx=zeros(length(ind),1);
% for i=1:length(ind)-1
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% xx(end)=0.5;
% figure
% plot(xx,zz)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',13);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',13);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',12.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
%------------------------------------------------------------------
%------------------------------------------------------------------
% %Comandi per plottare soluzione relativa a screenheiko1_1
% % ind = [28 26 25 31 12 10 9 15];
% ind = [28 25 12 9];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% xx=zeros(length(ind),1);
% for i=1:length(ind)
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% figure
% plot(xx,zz)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',13);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',13);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',12.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
%------------------------------------------------------------------
%------------------------------------------------------------------
% %Comandi per plottare soluzione relativa a screenheiko1_2
% %ind = [112 110 109 104 100 98 97 123 48 46 45 40 36 34 33 59];
% ind = [112 109 100 97 48 45 36 33];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% xx=zeros(length(ind),1);
% for i=1:length(ind)
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% figure
% plot(xx,zz)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',13);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',13);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',12.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
% figure
% l=find(zz<0,1);
% xxx = xx(1:l-1);
% xxx=abs(xxx-(-0.5));
% yyy=1./sqrt(xxx);
% loglog(xxx,zz(1:l-1),'-*')
% hold on
% loglog(xxx,yyy,'--')
% box off
% legend(strcat('$\phi_3$'),'$r^{-1/2}$','interpreter','latex','fontsize',13)
% xlabel('$r=|x_1-(-0.5)|$','interpreter','latex','fontsize',13);
% name=strcat('Comportamento_asintotico_phi3_sott_temporale_',num2str(hk));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
%------------------------------------------------------------------
%------------------------------------------------------------------
% %Comandi per plottare soluzione relativa a screenheiko1_3
% % ind = [448 446 445 440 436 434 433 416 400 398 397 392 388 386 385 491 192 190 ...
% %        189 184 180 178 177 160 144 142 141 136 132 130 129 235];
%    ind = [448 445 436 433 400 397 388 385 192 ...
%        189 180 177 144 141 132 129];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% xx=zeros(length(ind),1);
% for i=1:length(ind)
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% figure
% plot(xx,zz)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',13);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',13);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',12.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
% 
% figure
% l=find(zz<0,1);
% xxx = xx(1:l-1);
% xxx=abs(xxx-(-0.5));
% yyy=1./sqrt(xxx);
% loglog(xxx,zz(1:l-1),'-*')
% hold on
% loglog(xxx,yyy,'--')
% box off
% legend(strcat('$\phi_3$'),'$r^{-1/2}$','interpreter','latex','fontsize',13)
% xlabel('$r=|x_1-(-0.5)|$','interpreter','latex','fontsize',13);
% name=strcat('Comportamento_asintotico_phi3_sott_temporale_',num2str(hk));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
%------------------------------------------------------------------
%------------------------------------------------------------------
% %Comandi per plottare soluzione relativa a screenheiko1_4
% % ind = [1792 1790 1789 1784 1780 1778 1777 1760 1744 1742 1741 1736 1732 1730 1729 1664 ... 
% %        1600 1598 1597 1592 1588 1586 1585 1568 1552 1550 1549 1544 1540 1538 1537 1963 ...
% %         768 766 765 760 756 754 753 736 720 718 717 712 708 706 705 640 576 574 573 ...
% %          568 564 562 561 544 528 526 525 520 516 514 513 939];
% % ind = [1279 1271 1270 1269 1267 1247 1246 1245 1231 1223 1222 1221 1219 1151 1150 1149 ...
% %        1087 1079 1078 1077 1075 1055 1054 1053 1039 1031 1030 1029 1027 1449 1450 1452 ...
% %         255 247 246 245 243 223 222 221 207 199 198 197 195 127 126 125 63 55 54 53 51 ...
% %          31 30 29 15 7 6 5 3 425 426 428];
% ind = [1792 1789 1780 1777 1744 1741 1732 1729 ... 
%        1600 1597 1588 1585 1552 1549 1540 1537 ...
%         768 765 756 753 720 717 708 705 576 573 ...
%          564 561 528 525 516 513];
% 
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% xx=zeros(length(ind),1);
% for i=1:length(ind)
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% figure
% plot(xx,zz)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',13);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',13);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',12.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
% 
% figure
% l=find(zz<0,1);
% xxx = xx(1:l-1);
% xxx=abs(xxx-(-0.5));
% yyy=1./sqrt(xxx);
% loglog(xxx,zz(1:l-1),'-*')
% hold on
% loglog(xxx,yyy,'--')
% box off
% legend(strcat('$\phi_3$'),'$r^{-1/2}$','interpreter','latex','fontsize',13)
% xlabel('$r=|x_1-(-0.5)|$','interpreter','latex','fontsize',13);
% name=strcat('Comportamento_asintotico_phi3_sott_temporale_',num2str(hk));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
%------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% screengraded
%----------------------------------------------------------
% % %Screengraded 0
% % ind=[9 10 25 26 41 42 57 58 73 74 89 90 105 106 121 122];
% ind=[9 25 41 57 73 89 105 121];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% xx=zeros(length(ind),1);
% for i=1:length(ind)
%     cS(i,:) = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% figure
% plot(xx,zz,'linewidth',1.5)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',14);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',14);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',13.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
% 
% figure
% l=find(zz<0,1);
% xxx = xx(1:l-1);
% xxx=abs(xxx-(-0.5));
% yyy=1./sqrt(xxx);
% loglog(xxx,zz(1:l-1),'-*')
% hold on
% loglog(xxx,yyy,'--')
% box off
% legend(strcat('$\phi_3$'),'$r^{-1/2}$','interpreter','latex','fontsize',13)
% xlabel('$r=|x_1-(-0.5)|$','interpreter','latex','fontsize',13);
% name=strcat('Comportamento_asintotico_phi3_sott_temporale_',num2str(hk));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
%------------------------------------------------------------------
%------------------------------------------------------------------
% % %Screengraded 1
% % ind=[17 18 49 50 81 82 113 114 145 146 177 178 209 210 241 242 273 274 305 306 ...
% %     337 338 369 370 401 402 433 434 465 466 497 498];
% % ind=[3 4 35 36 67 68 99 100 131 132 163 164 195 196 227 228 259 260 291 292 323 324 ...
% %      355 356 387 388 419 420 451 452 483 484];
% ind=[17 49 81 113 145 177 209 241 273 305 337 369 401 433 465 497];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% xx=zeros(length(ind),1);
% for i=1:length(ind)
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% figure
% plot(xx,zz,'linewidth',1.5)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',14);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',14);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',13.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
%   
% figure
% l=find(zz<0,1);
% xxx = xx(1:l-1);
% xxx=abs(xxx-(-0.5));
% yyy=1./sqrt(xxx);
% loglog(xxx,zz(1:l-1),'-*')
% hold on
% loglog(xxx,yyy,'--')
% box off
% legend(strcat('$\phi_3$'),'$r^{-1/2}$','interpreter','latex','fontsize',13)
% xlabel('$r=|x_1-(-0.5)|$','interpreter','latex','fontsize',13);
% name=strcat('Comportamento_asintotico_phi3_sott_temporale_',num2str(hk));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
%------------------------------------------------------------------
%------------------------------------------------------------------
% %Screengraded 2
% ind=[33 97 161 225 289 353 417 481 545 609 673 737 801 865 929 993 1057 1121 1185 1249 ...
%      1313 1377 1441 1505 1569 1633 1697 1761 1825 1889 1953 2017];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% xx=zeros(length(ind),1);
% for i=1:length(ind)
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% figure
% plot(xx,zz,'linewidth',1.5)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',14);
% ylabel(strcat('$\phi_3(x_1,0,0,t)$'),'interpreter','latex','fontsize',14);
% title(strcat('$\phi_3(x_1,0,0,t)$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',13.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
% 
% figure
% l=find(zz<0,1);
% xxx = xx(1:l-1);
% xxx=abs(xxx-(-0.5));
% yyy=1./sqrt(xxx);
% loglog(xxx,yyy,'--','linewidth',1)
% hold on
% loglog(xxx,zz(1:l-1),'-*','linewidth',1)
% box off
% legend('$r^{(-1/2)}$',strcat('$\phi_3$'),'interpreter','latex','fontsize',13)
% xlabel('$r=|x_1-(-0.5)|$','interpreter','latex','fontsize',14);
% title(strcat('Comportamento asintotico di $\phi_3(x_1,0,0,',num2str(pb_param.T_fin),')$ per $r \to 0$'),'interpreter','latex','fontsize',13)
% name=strcat('Comportamento_asintotico_phi3_sott_temporale_',num2str(hk));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% screen1
%----------------------------------------------------------
% %screen1_0
% ind=[8 7 3 4];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% %xx = linspace(-0.5,0.5,length(ind));
% xx=zeros(length(ind),1);
% for i=2:length(ind)-1
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% xx(1)=-0.5;
% xx(end)=0.5;
% figure
% plot(xx,zz)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',13);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',13);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',12.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
%-----------------------------------------------------------
% %screen1_1
% %ind=[32 28 26 25 9 10 11 13];
% ind=[32 28 25 9 11 13];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% %xx = linspace(-0.5,0.5,length(ind));
% xx=zeros(length(ind),1);
% for i=2:length(ind)-1
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% xx(1)=-0.5;
% xx(end)=0.5;
% figure
% plot(xx,zz,'linewidth',1.5)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',14);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',14);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',13.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
%-----------------------------------------------------------
%-----------------------------------------------------------
% %screen1_2
% ind=[128 112 109 100 97 33 35 41 43 49];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% %xx = linspace(-0.5,0.5,length(ind));
% xx=zeros(length(ind),1);
% for i=2:length(ind)-1
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% xx(1)=-0.5;
% xx(end)=0.5;
% figure
% plot(xx,zz,'linewidth',1.5)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',14);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',14);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',13.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
%------------------------------------------------------------------
%------------------------------------------------------------------
% %screen1_3
% % ind=[512 448 446 445 440 436 434 433 416 400 398 397 392 388 386 385 ...
% %      129 130 131 133 137 138 139 145 161 162 163 165 169 170 171 193];
% ind=[512 448 445 436 433 400 397 388 385 129 131 137 139 161  163 169 171 193];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% xx=zeros(length(ind),1);
% for i=2:length(ind)-1
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% xx(1)=-0.5;
% xx(end)=0.5;
% figure
% plot(xx,zz,'linewidth',1.5)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',14);
% ylabel(strcat('$\phi_3$'),'interpreter','latex','fontsize',14);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',13.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
%------------------------------------------------------------------
%------------------------------------------------------------------
% % screen1_4
% %  ind=[2048 1792 1790 1789 1784 1780 1778 1777 1760 1744 1742 1741 1736 1732 ...
% %      1730 1729 1664 1600 1598 1597 1592 1588 1586 1585 1568 1552 1550 1549 1544 1540 ...
% %       1538 1537 513 514 515 517 521 522 523 529 545 546 547 549 553 554 555 577 641 ...
% %        642 643 645 649 650 651 657 673 674 675 677 681 682 683 769];
% % ind=[1280 1536 1534 1535 1527 1532 1530 1531 1499 1520 1518 1519 1511 1516 1514 1515 1387 ...
% %     1472 1470 1471 1463 1468 1466 1467 1435 1456 1454 1455 1447 1452 1450 1451 427 426 425 ...
% %     421 419 418 417 401 395 394 393 389 387 386 385 321 299 298 297 293 291 290 289 273 267 ...
% %     266 265 261 259 258 257 171];
% % ind=[1028 1032  1030 1031 1036 1044 1042 1043 1060 1064 1062 1063 1068 1092 1090 1091 1156 ...
% %      1160 1158 1159 1164 1172 1170 1171 1188 1192 1190 1191 1196 1284 1282 1283 511 510 ...
% %       509 255 247 246 245 243 223 221 207 199 198 197 195 127 126 125 63 55 54 53 51 31 ...
% %       30 29 15 7 6 5 3];
% ind=[2048 1792 1789 1780 1777 1744 1741 1732 1729 1600 1597 1588 1585 1552 1549 1540 ...
%       1537 513 515 521 523 545 547 553 555 641 643 649 651 673 675 681 683 769];
% hk=pb_param.Nt;
% zz = density(3*ind,hk);
% xx=zeros(length(ind),1);
% for i=2:length(ind)-1
%     cS = domainMesh.center(ind(i),:);
%     xx(i)=cS(1);
% end
% xx(1)=-0.5;
% xx(end)=0.5;
% figure
% plot(xx,zz,'linewidth',1.5)
% box off
% grid on
% xlabel('Asse $x_1$','interpreter','latex','fontsize',14);
% ylabel('$\phi_3$','interpreter','latex','fontsize',14);
% title(strcat('$\phi_3$ per $x_1\in[-0.5,0.5]$ e $ t =',num2str(pb_param.T_fin),'$'),'interpreter','latex','fontsize',13.5); 
% name=strcat(pb_param.domain_type,'_',num2str(pb_param.lev),'_Bicchierino_phi3','_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')
% 
% figure
% l=find(zz<0,1);
% xxx = xx(1:l-1);
% xxx=abs(xxx-(-0.5));
% yyy=1./sqrt(xxx);
% loglog(xxx,zz(1:l-1),'-*')
% hold on
% loglog(xxx,yyy,'--')
% box off
% legend(strcat('$\phi_3$'),'$r^{-1/2}$','interpreter','latex','fontsize',13,'location','best')
% xlabel('$r=|x_1-(-0.5)|$','interpreter','latex','fontsize',13);
% name=strcat('Comportamento_asintotico_phi3_t_',num2str(pb_param.T_fin));
% print(name,'-dtiff')
% savefig(strcat(name,'.fig'))
% print(name,'-depsc2','-r500')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Densità al variare del tempo in un baricentro di un triangolo 
% T=pb_param.T_fin;
% ind=1;
% for i=1:3
%     alpha = density(3*ind-(3-i),:);
%     tt = linspace(0,pb_param.T_fin,pb_param.Nt+1);
%     XX = linspace(0,pb_param.T_fin,2001);
%     YY=zeros(1,length(XX));
%     for l=1:pb_param.Nt-1
%         chi_l=@(t)((t >= tt(l) & t<tt(l+1)).*1);
%         YY=YY+chi_l(XX)*alpha(l);
%     end
%     chi_l=@(t)((t >= tt(pb_param.Nt) & t<=tt(pb_param.Nt+1)).*1);
%     YY=YY+chi_l(XX)*alpha(pb_param.Nt);
%     figure(i)
%     box off
%     plot(XX,YY)
%     xlabel('Asse $t$','interpreter','latex','fontsize',13);
%     ylabel(strcat('$\phi_',num2str(i)','(x_1,x_2,x_3,t)$'),'interpreter','latex','fontsize',16);
%     title(strcat('Grafico di $\phi_',num2str(i),'(x_1,x_2,x_3,t)$ per $t\in[0,',num2str(T),']$ sul triangolo $',num2str(ind),'$'),'interpreter','latex','fontsize',12); 
%     name=strcat('phi',num2str(i),'_punto_triangolo_',num2str(ind));
%     print(name,'-dtiff')
%     savefig(strcat(namae,'.fig'))
%     close all
% end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Post processing
% x=[0,0,2.75]; %Punto in cui si vuole valutare la soluzione u
% T=pb_param.T_fin;
% t=linspace(0,T,697);
% U=zeros(3,length(t));
% %Richiamiamo la funzione time3D_history_post_processing() per ricostruire
% %la soluzione u attraverso la fase di post-processing  
% for i=1:length(t)
%     U(:,i)=time3D_history_post_processing(pb_param,domainMesh,density,x,t(i));
% end
% %Grafico delle componenti di u al variare del tempo
% for j=1:3
%     figure(j)
%     plot(t,U(j,:))
%     box off
%     xlabel('Asse $t$','interpreter','latex','fontsize',13);
%     ylabel(strcat('$u_',num2str(j)','$'),'interpreter','latex','fontsize',16);
%     title(strcat('Grafico di $u_',num2str(j),'(',num2str(x(1)),',',num2str(x(2)),',',num2str(x(3)),',t)$ per $t\in[0,',num2str(T),']$'),'interpreter','latex','fontsize',13); 
%     xlim([0 36])
%     name=strrep(strcat('u',num2str(j),'_punto_','(',num2str(x(1)),';',num2str(x(2)),';',num2str(x(3)),')'),'.',',');
%     %print(name,'-dtiff')
%     savefig(strcat(name,'.fig'))
% end
% name=strrep(strcat('sol_punto_','(',num2str(x(1)),';',num2str(x(2)),';',num2str(x(3)),')'),'.',',');
% save(name)


% EL2_0=sqrt(trapz(ttt,(abs(U0(3,:)-sol_bar).^2)));
% EL2_1=sqrt(trapz(ttt,(abs(U1(3,:)-sol_bar).^2)));
% EL2_2=sqrt(trapz(ttt,(abs(U2(3,:)-sol_bar).^2)));
% EL2_3=sqrt(trapz(ttt,(abs(U3(3,:)-sol_bar).^2)));
% loglog([1 0.5 0.25 0.125],[1 0.5 0.25 0.125].^(3/2));
% hold on
% loglog([1 0.5 0.25 0.125],[EL2_0 EL2_1 EL2_2 EL2_3]);
return


