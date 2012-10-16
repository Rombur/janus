/*
Copyright (c) 2012, Bruno Turcksin.

This file is part of Janus.

Janu is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
he Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Janus is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Janus.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "LS.hh"

LS::LS(unsigned int sn,unsigned int L_max,bool galerkin) :
  QUADRATURE(sn,L_max,galerkin)
{}

void LS::Build_octant()
{
  switch (sn)
  {
    case 2:
      {
        const double w(1.);
        const double direction(0.577350269189625764509149);

        omega[0](0) = direction;
        omega[0](1) = direction;
        omega[0](2) = direction;

        weight(0) = w;

        break;
      }
    case 4:
      {
        const double w(1./3.);
        d_vector direction(2,0.);
        direction[0] = 0.350021174581540677777041;
        direction[1] = 0.868890300722201205229788;

        Compute_omega(direction);

        weight(0) = w;
        weight(1) = w;
        weight(2) = w;

        break;
      }
    case 6 :
      {
        const double w_1(0.176126130863383433783565);
        const double w_2(0.157207202469949899549768);
        d_vector direction(3,0.);
        
        direction[0] = 0.266635401516704720331535;
        direction[1] = 0.681507726536546927403750;
        direction[2] = 0.926180935517489107558380;

        Compute_omega(direction);

        weight(0) = w_1;
        weight(1) = w_2;
        weight(2) = w_1;
        weight(3) = w_2;
        weight(4) = w_2;
        weight(5) = w_1;

        break;
      }
    case 8:
      {
        const double w_1(0.120987654320987654320988);
        const double w_2(0.0907407407407407407407407);
        const double w_3(0.0925925925925925925925926);
        d_vector direction(4,0.);

        direction[0] = 0.218217890235992381266097;
        direction[1] = 0.577350269189625764509149;
        direction[2] = 0.786795792469443145800830;
        direction[3] = 0.951189731211341853132399;

        Compute_omega(direction);

        weight(0) = w_1;
        weight(1) = w_2;
        weight(2) = w_2;
        weight(3) = w_1;
        weight(4) = w_2;
        weight(5) = w_3;
        weight(6) = w_2;
        weight(7) = w_2;
        weight(8) = w_2;
        weight(9) = w_1;

        break;
      }
    case 10:
      {
        const double w_1(0.0893031479843567214704325);
        const double w_2(0.0725291517123655242296233);
        const double w_3(0.0450437674364086390490892);
        const double w_4(0.0539281144878369243545650);
        d_vector direction(5,0.);

        direction[0] = 0.189321326478010476671494;
        direction[1] = 0.508881755582618974382711;
        direction[2] = 0.694318887594384317279217;
        direction[3] = 0.839759962236684758403029;
        direction[4] = 0.963490981110468484701598;

        Compute_omega(direction);

        weight(0) = w_1; 
        weight(1) = w_2;
        weight(2) = w_3;
        weight(3) = w_2;
        weight(4) = w_1;
        weight(5) = w_2;
        weight(6) = w_4;
        weight(7) = w_4;
        weight(8) = w_2;
        weight(9) = w_3;
        weight(10) = w_4;
        weight(11) = w_3;
        weight(12) = w_2;
        weight(13) = w_2;
        weight(14) = w_1;

        break;
      }
    case 12:
      {
        const double w_1(0.0707625899700910439766549);
        const double w_2(0.0558811015648888075828962);
        const double w_3(0.0373376737588285824652402);
        const double w_4(0.0502819010600571181385765);
        const double w_5(0.0258512916557503911218290);
        d_vector direction(6,0.);

        direction[0] = 0.167212652822713264084504; 
        direction[1] = 0.459547634642594690016761; 
        direction[2] = 0.628019096642130901034766; 
        direction[3] = 0.760021014833664062877138; 
        direction[4] = 0.872270543025721502340662; 
        direction[5] = 0.971637719251358378302376; 

        Compute_omega(direction);

        weight(0) = w_1;
        weight(1) = w_2;
        weight(2) = w_3;
        weight(3) = w_3;
        weight(4) = w_2;
        weight(5) = w_1;
        weight(6) = w_2;
        weight(7) = w_4;
        weight(8) = w_5;
        weight(9) = w_4;
        weight(10) = w_2;
        weight(11) = w_3;
        weight(12) = w_5;
        weight(13) = w_5;
        weight(14) = w_3;
        weight(15) = w_3;
        weight(16) = w_4;
        weight(17) = w_3;
        weight(18) = w_2;
        weight(19) = w_2;
        weight(20) = w_1;

        break;
      }
    case 14:
      {
        const double w_1(1.0579970408969969964063611); 
        const double w_2(0.0489007976368104874582568); 
        const double w_3(0.0227935342411872473257345); 
        const double w_4(0.0394132005950078294492985); 
        const double w_5(0.0380990861440121712365891); 
        const double w_6(0.0258394076418900119611012); 
        const double w_7(0.00826957997262252825269908);
        d_vector direction(7,0.);
        
        direction[0] = 0.151985861461031912404799;
        direction[1] = 0.422156982304796966896263;
        direction[2] = 0.577350269189625764509149;
        direction[3] = 0.698892086775901338963210;
        direction[4] = 0.802226255231412057244328;
        direction[5] = 0.893691098874356784901111;
        direction[6] = 0.976627152925770351762946;

        Compute_omega(direction);

        weight(0) = w_1;
        weight(1) = w_2;
        weight(2) = w_3;
        weight(3) = w_4;
        weight(4) = w_3;
        weight(5) = w_2;
        weight(6) = w_1;
        weight(7) = w_2;
        weight(8) = w_5;
        weight(9) = w_6;
        weight(10) = w_6;
        weight(11) = w_5;
        weight(12) = w_2;
        weight(13) = w_3;
        weight(14) = w_6;
        weight(15) = w_7;
        weight(16) = w_6;
        weight(17) = w_3;
        weight(18) = w_4;
        weight(19) = w_6;
        weight(20) = w_6;
        weight(21) = w_4;
        weight(22) = w_3;
        weight(23) = w_5;
        weight(24) = w_3;
        weight(25) = w_2;
        weight(26) = w_2;
        weight(27) = w_1;
        break;
      }
    case 16:
      {
        const double w_1(0.0489872391580385335008367);
        const double w_2(0.0413295978698440232405505);
        const double w_3(0.0203032007393652080748070);
        const double w_4(0.0265500757813498446015484);
        const double w_5(0.0379074407956004002099321);
        const double w_6(0.0135295047786756344371600);
        const double w_7(0.0326369372026850701318409);
        const double w_8(0.0103769578385399087825920);
        d_vector direction(8,0.);

        direction[0] = 0.138956875067780344591732;
        direction[1] = 0.392289261444811712294197;
        direction[2] = 0.537096561300879079878296;
        direction[3] = 0.650426450628771770509703;
        direction[4] = 0.746750573614681064580018;
        direction[5] = 0.831996556910044145168291;
        direction[6] = 0.909285500943725291652116;
        direction[7] = 0.980500879011739882135849;

        Compute_omega(direction);

        weight(0) = w_1; 
        weight(1) = w_2;
        weight(2) = w_3;
        weight(3) = w_4;
        weight(4) = w_4;
        weight(5) = w_3;
        weight(6) = w_2;
        weight(7) = w_1;
        weight(8) = w_2;
        weight(9) = w_5;
        weight(10) = w_6;
        weight(11) = w_7;
        weight(12) = w_6;
        weight(13) = w_5;
        weight(14) = w_2;
        weight(15) = w_3;
        weight(16) = w_6;
        weight(17) = w_8;
        weight(18) = w_8;
        weight(19) = w_6;
        weight(20) = w_3;
        weight(21) = w_4;
        weight(22) = w_7;
        weight(23) = w_8;
        weight(24) = w_7;
        weight(25) = w_4;
        weight(26) = w_4;
        weight(27) = w_6;
        weight(28) = w_6;
        weight(29) = w_4;
        weight(30) = w_3;
        weight(31) = w_5;
        weight(32) = w_3;
        weight(33) = w_2;
        weight(34) = w_2;
        weight(35) = w_1;
        
        break;
      }
    case 18:
      {
        const double w_1(0.0422646448843821748535825);  
        const double w_2(0.0376127473827281471532380);  
        const double w_3(0.0122691351637405931037187);  
        const double w_4(0.0324188352558815048715646);  
        const double w_5(0.00664438614619073823264082); 
        const double w_6(0.0312093838436551370068864);  
        const double w_7(0.0160127252691940275641645);  
        const double w_8(0.0200484595308572875885066);  
        const double w_9(0.000111409402059638628382279);
        const double w_10(0.0163797038522425240494567);
        d_vector direction(9,0.);

        direction[0] = 0.129344504545924818514086;
        direction[1] = 0.368043816053393605686086;
        direction[2] = 0.504165151725164054411848;
        direction[3] = 0.610662549934881101060239;
        direction[4] = 0.701166884252161909657019;
        direction[5] = 0.781256199495913171286914;
        direction[6] = 0.853866206691488372341858;
        direction[7] = 0.920768021061018932899055;
        direction[8] = 0.983127661236087115272518;

        Compute_omega(direction);
        
        weight(0) = w_1;
        weight(1) = w_2;
        weight(2) = w_3;
        weight(3) = w_4;
        weight(4) = w_5;
        weight(5) = w_4;
        weight(6) = w_3;
        weight(7) = w_2;
        weight(8) = w_1;
        weight(9) = w_2;
        weight(10) = w_6;
        weight(11) = w_7;
        weight(12) = w_8;
        weight(13) = w_8;
        weight(14) = w_7;
        weight(15) = w_6;
        weight(16) = w_2;
        weight(17) = w_7;
        weight(18) = w_9;
        weight(19) = w_10;
        weight(20) = w_9;
        weight(21) = w_7;
        weight(22) = w_3;
        weight(23) = w_4;
        weight(24) = w_8;
        weight(25) = w_10;
        weight(26) = w_10;
        weight(27) = w_8;
        weight(28) = w_4;
        weight(29) = w_5;
        weight(30) = w_8;
        weight(31) = w_9;
        weight(32) = w_8;
        weight(33) = w_5;
        weight(34) = w_4;
        weight(35) = w_7;
        weight(36) = w_7;
        weight(37) = w_4;
        weight(38) = w_3;
        weight(39) = w_6;
        weight(40) = w_3;
        weight(41) = w_2;
        weight(42) = w_2;
        weight(43) = w_1;

        break;
      }
    case 20:
      {
        const double w_1(0.0370210490604481342320295);   
        const double w_2(0.0332842165376314841003910);   
        const double w_3(0.0111738965965092519614021);   
        const double w_4(0.0245177476959359285418987);   
        const double w_5(0.0135924329650041789567081);   
        const double w_6(0.0318029065936585971501960);   
        const double w_7(0.00685492401402507781062634);  
        const double w_8(0.0308105481755299327227893);
        const double w_9(-0.000139484716502602877593527);
        const double w_10(0.00544675187330776223879437);
        const double w_11(0.00474564692642379971238396);
        const double w_12(0.0277298541009064049325246);
        d_vector direction(10,0.);

        direction[0] = 0.120603343036693597409418;
        direction[1] = 0.347574292315847257336779;
        direction[2] = 0.476519266143665680817278;
        direction[3] = 0.577350269489625764509149;
        direction[4] = 0.663020403653288019308789;
        direction[5] = 0.738822561910371432904974;
        direction[6] = 0.807540401661143067193530;
        direction[7] = 0.870852583760463975580977;
        direction[8] = 0.929863938955324566667817;
        direction[9] = 0.985347485558646574628509;

        Compute_omega(direction);

        weight(0) = w_1;
        weight(1) = w_2;
        weight(2) = w_3;
        weight(3) = w_4;
        weight(4) = w_5;
        weight(5) = w_5;
        weight(6) = w_4;
        weight(7) = w_3;
        weight(8) = w_3;
        weight(9) = w_1;
        weight(10) = w_2;
        weight(11) = w_6;
        weight(12) = w_7;
        weight(13) = w_8;
        weight(14) = w_9;
        weight(15) = w_8;
        weight(16) = w_7;
        weight(17) = w_6;
        weight(18) = w_2;
        weight(19) = w_3;
        weight(20) = w_7;
        weight(21) = w_10;
        weight(22) = w_11;
        weight(23) = w_11;
        weight(24) = w_10;
        weight(25) = w_7;
        weight(26) = w_3;
        weight(27) = w_4;
        weight(28) = w_8;
        weight(29) = w_11;
        weight(30) = w_12;
        weight(31) = w_11;
        weight(32) = w_8;
        weight(33) = w_4;
        weight(34) = w_5;
        weight(35) = w_9;
        weight(36) = w_11;
        weight(37) = w_11;
        weight(38) = w_9;
        weight(39) = w_5;
        weight(40) = w_5;
        weight(41) = w_8;
        weight(42) = w_10;
        weight(43) = w_8;
        weight(44) = w_5;
        weight(45) = w_4;
        weight(46) = w_7;
        weight(47) = w_7;
        weight(48) = w_4;
        weight(49) = w_3;
        weight(50) = w_6;
        weight(51) = w_3;
        weight(52) = w_2;
        weight(53) = w_2;
        weight(54) = w_1;

        break;
      }
    case 22:
      {
        const double w_1(0.0329277718552552308051381);   
        const double w_2(0.0309569328165031538543025);   
        const double w_3(0.00577105953220643022391829);  
        const double w_4(0.0316834548379952775919418);
        const double w_5(0.00669350304140992494103696);
        const double w_6(0.0368381622687682466526634);   
        const double w_7(0.0273139698006629537455404);   
        const double w_8(0.0100962716435030437817055);   
        const double w_9(0.0195181067555849392224199);   
        const double w_10(0.0117224275470949786864925);  
        const double w_11(-0.00442773155233893239996431);
        const double w_12(0.0156214785078803432781324);
        const double w_13(-0.0101774221315738297143270); 
        const double w_14(0.0135061258938431808485310);
        d_vector direction(11,0.);

        direction[0] = 0.113888641383070838173488;
        direction[1] = 0.330271760593086736334651;
        direction[2] = 0.452977095507524183904005;
        direction[3] = 0.548905330875560154226714;
        direction[4] = 0.630401360620980621392149;
        direction[5] = 0.702506006153654989703184;
        direction[6] = 0.767869456282208576047898;
        direction[7] = 0.828089557415325768804621;
        direction[8] = 0.884217805921983001958912;
        direction[9] = 0.936989829997455780115072;
        direction[10] = 0.986944149751056870330152;

        Compute_omega(direction);

        weight(0) = w_1;
        weight(1) = w_2;
        weight(2) = w_3;
        weight(3) = w_4;
        weight(4) = w_5;
        weight(5) = w_6;
        weight(6) = w_5;
        weight(7) = w_4;
        weight(8) = w_3;
        weight(9) = w_3;
        weight(10) = w_1;
        weight(11) = w_2;
        weight(12) = w_7;
        weight(13) = w_8;
        weight(14) = w_9;
        weight(15) = w_10;
        weight(16) = w_10;
        weight(17) = w_9;
        weight(18) = w_8;
        weight(19) = w_7;
        weight(20) = w_2;
        weight(21) = w_3;
        weight(22) = w_8;
        weight(23) = w_11;
        weight(24) = w_12;
        weight(25) = w_13;
        weight(26) = w_12;
        weight(27) = w_11;
        weight(28) = w_8;
        weight(29) = w_3;
        weight(30) = w_4;
        weight(31) = w_9;
        weight(32) = w_12;
        weight(33) = w_14;
        weight(34) = w_14;
        weight(35) = w_12;
        weight(36) = w_9;
        weight(37) = w_4;
        weight(38) = w_5;
        weight(39) = w_10;
        weight(40) = w_13;
        weight(41) = w_14;
        weight(42) = w_13;
        weight(43) = w_10;
        weight(44) = w_5;
        weight(45) = w_6;
        weight(46) = w_10;
        weight(47) = w_12;
        weight(48) = w_12;
        weight(49) = w_10;
        weight(50) = w_6;
        weight(51) = w_5;
        weight(52) = w_9;
        weight(53) = w_11;
        weight(54) = w_9;
        weight(55) = w_5;
        weight(56) = w_4;
        weight(57) = w_8;
        weight(58) = w_8;
        weight(59) = w_4;
        weight(60) = w_3;
        weight(61) = w_7;
        weight(62) = w_3;
        weight(63) = w_2;
        weight(64) = w_2;
        weight(65) = w_1;

        break;
      }
    case 24:
      {
        const double w_1(0.0295284942030736546025272);   
        const double w_2(0.0281530651743695026834932);   
        const double w_3(0.00519730128072174996473824);  
        const double w_4(0.0259897467786242920448933);   
        const double w_5(0.00146378160153344429844948);  
        const double w_6(0.0166609651269037212368055);   
        const double w_7(0.0281343344093849194875108);   
        const double w_8(0.00214364311909247909952968);  
        const double w_9(0.0331943417648083019611294);   
        const double w_10(-0.0142483904822400753741381); 
        const double w_11(0.0416812529998231580614934);  
        const double w_12(0.00323522898964475022578598); 
        const double w_13(0.000813552611571786631179287);
        const double w_14(0.00228403610697848813660369); 
        const double w_15(0.0338971925236628645848122);  
        const double w_16(-0.00644725595698339499416262);
        d_vector direction(12,0.);

        direction[0] = 0.107544208775147285552086;
        direction[1] = 0.315151630853896874875332;
        direction[2] = 0.432522073446742487657060;
        direction[3] = 0.524242441631224399254880;
        direction[4] = 0.602150256328323868809286;
        direction[5] = 0.671073561381361944701265;
        direction[6] = 0.733549261041044861004094;
        direction[7] = 0.791106384731321324814121; 
        direction[8] = 0.844750913317919895113069;
        direction[9] = 0.895186516397704814461305;
        direction[10] = 0.942928254285052510917188;
        direction[11] = 0.988366574868785749937406;

        Compute_omega(direction);

        weight(0) = w_1;
        weight(1) = w_2;
        weight(2) = w_3;
        weight(3) = w_4;
        weight(4) = w_5;
        weight(5) = w_6;
        weight(6) = w_6;
        weight(7) = w_5;
        weight(8) = w_4;
        weight(9) = w_3;
        weight(10) = w_2;
        weight(11) = w_1;
        weight(12) = w_2;
        weight(13) = w_7;
        weight(14) = w_8;
        weight(15) = w_9;
        weight(16) = w_10;
        weight(17) = w_11;
        weight(18) = w_10;
        weight(19) = w_9;
        weight(20) = w_8;
        weight(21) = w_7;
        weight(22) = w_2;
        weight(23) = w_3;
        weight(24) = w_8;
        weight(25) = w_12;
        weight(26) = w_13;
        weight(27) = w_14;
        weight(28) = w_14;
        weight(29) = w_13;
        weight(30) = w_12;
        weight(31) = w_8;
        weight(32) = w_3;
        weight(33) = w_4;
        weight(34) = w_9;
        weight(35) = w_13;
        weight(36) = w_15;
        weight(37) = w_16;
        weight(38) = w_15;
        weight(39) = w_13;
        weight(40) = w_9;
        weight(41) = w_4;
        weight(42) = w_5;
        weight(43) = w_10;
        weight(44) = w_14;
        weight(45) = w_16;
        weight(46) = w_16;
        weight(47) = w_14;
        weight(48) = w_10;
        weight(49) = w_5;
        weight(50) = w_6;
        weight(51) = w_11;
        weight(52) = w_14;
        weight(53) = w_15;
        weight(54) = w_14;
        weight(55) = w_11;
        weight(56) = w_6;
        weight(57) = w_6;
        weight(58) = w_10;
        weight(59) = w_13;
        weight(60) = w_13;
        weight(61) = w_10;
        weight(62) = w_6;
        weight(63) = w_5;
        weight(64) = w_9;
        weight(65) = w_12;
        weight(66) = w_9;
        weight(67) = w_5;
        weight(68) = w_4;
        weight(69) = w_8;
        weight(70) = w_8;
        weight(71) = w_4;
        weight(72) = w_3;
        weight(73) = w_7;
        weight(74) = w_3;
        weight(75) = w_2;
        weight(76) = w_2;
        weight(77) = w_1;

        break;
      }
    default :
      {
        Check(false,"This quadrature ordr does not exist.");
      }
  } 
}

void LS::Compute_omega(d_vector const &direction)
{
  unsigned int offset(0);
  const unsigned int size(direction.size());

  for (unsigned int j=0; j<size; ++j)
  {
    for (unsigned int i=0; i<size-j; ++i)
    {
      omega[offset+i](0) = direction[size-1-j-i];
      omega[offset+i](1) = direction[i];
      omega[offset+i](2) = direction[j];
    }
    offset += size-j;
  }
}
