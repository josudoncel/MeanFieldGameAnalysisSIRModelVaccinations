function switchingcurve_mfopt_mfe_dynamics_ms04_mi04()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab code to generate the figure where we compare the switching
%   curve of the mean-field optimum and of the mean-field 
%   equilibrium. We also illustrate the dynamics for both
%   strategies when at time 0 the proportion of susceptible
%   and infected population are both 0.4.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=0.01;
addpath('ternaryplotfiles')

% switching curve of the mean-field optimum (obtained from global_optimum.py with different values of S0 and I0)
x=[0 0.02 0.04 0.06 0.08 0.1 0.12 0.145 0.17 0.195...
    0.22 0.245 0.27 0.3 0.33 0.365 0.4 0.435 0.47 0.655 1];
lin1=0.35:-0.02:0.01;
lin1=[lin1 0 0 0];
% switching curve of the mean-field equilibrium (obtained using equilibrium.py with different values of S0 and I0)
lin2_=0.35-0.54.*x-lin1;
lin2=max(0,lin2_);
lin3=-x+1-lin2-lin1;
% dynamics for the mean-field equilibrium (obtained using equilibrium.py with input S0=I0=0.4)
tray_eq=[0.4 0.4
0.38432 0.39708
0.369336577651 0.393726802349
0.355027719066 0.389971266872
0.341370547391 0.385844210116
0.328341584924 0.381376153439
0.315916988576 0.376597104338
0.304072758803 0.371536369916
0.292784922718 0.366222400911
0.282029692485 0.360682664284
0.271783600435 0.354943542163
0.26202361256 0.349030254744
0.252727222153 0.342966804728
0.243872525438 0.336775940848
0.235438281042 0.330479138149
0.227403955097 0.324096592742
0.219749753713 0.317647228939
0.212456644452 0.311148716807
0.20550636831 0.304617498341
0.198881443613 0.298068820666
0.19256516309 0.291516774798
0.186541585264 0.284974338713
0.180795521173 0.278453423588
0.175312517337 0.271964922252
0.170078835738 0.265518759015
0.165081431518 0.259123940173
0.161958743296 0.252788604579
0.158970026592 0.246550537216
0.156108853974 0.240412615226
0.153369124711 0.234377284033
0.150745049271 0.228446588606
0.148231134236 0.222622203156
0.145822167694 0.216905459284
0.143513205123 0.211297372591
0.141299555811 0.205798667803
0.139176769806 0.200409802433
0.137140625413 0.195130989037
0.135187117244 0.189962216106
0.133312444805 0.184903267658
0.131513001618 0.179953741575
0.129785364879 0.175113066747
0.12812628561 0.170380519079
0.126532679328 0.165755236415
0.125001617173 0.161236232441
0.123530317517 0.156822409612
0.122116138007 0.152512571171
0.120756568036 0.148305432296
0.119449221619 0.144199630434
0.118191830665 0.140193734876
0.11698223861 0.136286255608
0.115818394408 0.132475651481
0.114698346848 0.128760337762
0.113620239203 0.125138693079
0.112582304161 0.121609065823
0.111582859056 0.118169780026
0.110620301357 0.114819140754
0.109693104416 0.111555439057
0.108799813459 0.108376956488
0.107939041796 0.10528196924
0.10710946725 0.102268751909
0.106309828788 0.099335580926
0.10553892334 0.0964807376701
0.104795602798 0.093702511287
0.104078771184 0.0909992012393
0.103387381976 0.0883691196022
0.102720435585 0.0858105931271
0.102076976975 0.0833219650878
0.101456093411 0.0809015969268
0.100856912332 0.0785478697172
0.100278599353 0.076259185452
0.0997203563587 0.0740339681773
0.0991814197194 0.0718706649781
0.0986610585944 0.0697677468314
0.0981585733321 0.0677237093344
0.0976732939558 0.06573707332
0.0972045787323 0.0638063853673
0.0967518128171 0.0619302182165
0.0963144069729 0.0601071710959
0.0958917963555 0.0583358699683
0.095483439366 0.0566149677039
0.095088816562 0.0549431441868
0.0947074296272 0.0533191063587
0.0943388003947 0.0517415882091
0.0939824699213 0.0502093507129
0.0936379976094 0.0487211817238
0.0933049603748 0.0472758958254
0.092982951857 0.0458723341456
0.0926715816693 0.044509364137
0.0923704746876 0.0431858793277
0.0920792703749 0.0419007990449
0.0917976221396 0.0406530681151
0.0915251967257 0.0394416565428
0.0912616736334 0.0382655591713
0.0910067445684 0.0371237953266
0.090760112918 0.0360154084475
0.0905214932527 0.0349394657044
0.0902906108523 0.0338950576067
0.090067201254 0.0328812976023];
% dynamics for the mean-field optimum (obtained using global_optimum.py with input S0=I0=0.4)
tray_opt=[0.4 0.4
0.38432 0.39708
0.369336577651 0.393726802349
0.355027719066 0.389971266872
0.341370547391 0.385844210116
0.328341584924 0.381376153439
0.315916988576 0.376597104338
0.304072758803 0.371536369916
0.292784922718 0.366222400911
0.282029692485 0.360682664284
0.271783600435 0.354943542163
0.26202361256 0.349030254744
0.252727222153 0.342966804728
0.243872525438 0.336775940848
0.235438281042 0.330479138149
0.227403955097 0.324096592742
0.219749753713 0.317647228939
0.212456644452 0.311148716807
0.20550636831 0.304617498341
0.198881443613 0.298068820666
0.19256516309 0.291516774798
0.186541585264 0.284974338713
0.180795521173 0.278453423588
0.175312517337 0.271964922252
0.170078835738 0.265518759015
0.165081431518 0.259123940173
0.160307928981 0.252788604579
0.155746596401 0.246520073802
0.151386320079 0.240324901466
0.147216577999 0.234208921442
0.143227413395 0.228177294634
0.13940940847 0.22223455417
0.135753658481 0.216384648848
0.133609282917 0.210630984728
0.131554898315 0.204997338388
0.129586204822 0.19948362903
0.127699130994 0.194089550398
0.125889821109 0.188814591694
0.124154623141 0.183658057064
0.122490077371 0.178619083752
0.120892905597 0.173696658969
0.11936000095 0.168889635563
0.117888418255 0.164196746561
0.116475364939 0.159616618627
0.115118192454 0.155147784532
0.11381438818 0.150788694671
0.112561567799 0.146537727696
0.111357468104 0.14239320033
0.110199940227 0.138353376396
0.109086943259 0.134416475125
0.108016538244 0.130580678798
0.106986882523 0.126844139743
0.10599622441 0.123204986755
0.105042898181 0.119661330968
0.104125319351 0.116211271217
0.103241980243 0.112852898926
0.102391445799 0.109584302559
0.101572349652 0.106403571663
0.100783390414 0.103308800535
0.100023328198 0.100298091531
0.0992909813267 0.0973695580618
0.0985852232418 0.0945213272775
0.0979049795928 0.0917515424809
0.0972492254915 0.0890583652816
0.096616982927 0.0864399775132
0.0960073183293 0.0838945829317
0.0954193402725 0.0814204087116
0.0948521973095 0.0790157067565
0.0943050759308 0.0766787548387
0.0937771986376 0.0744078575803
0.0932678221254 0.0722013472908
0.092776235569 0.070057584671
0.0923017590036 0.067974959396
0.0918437417962 0.0659518905854
0.0914015612022 0.063986827173
0.0909746210015 0.0620782481819
0.0905623502092 0.0602246629156
0.0901642018572 0.0584246110712];

ternplot(tray_eq(1:size(tray_opt,1),1),tray_eq(1:size(tray_opt,1),2))
hold on
ternlabel('Proportion of Infected Population',...
    'Proportion of Susceptible Population',...
    'Proportion of Vaccinated Population')

for i=1:size(x,2)
    [lin1_x(i) lin1_y_(i)]=terncoords(x(i),lin1(i));
end
x=[0,lin1_x];
lin1_y=[0,lin1_y_];

k=1;
lin2_y_=lin1_y(2)-0.63*(x(2:size(x,2))-x(2));
lin2_y=[0 lin2_y_];
lin2_y=max(0,lin2_y);


area(x,[lin1_y;lin2_y-lin1_y]')

ternplot(tray_eq(1:size(tray_opt,1),1),tray_eq(1:size(tray_opt,1),2))
ternplot(tray_opt(1:size(tray_opt,1),1),tray_opt(1:size(tray_opt,1),2))


end
