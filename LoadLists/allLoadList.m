loadListName = 'allLoadList'; 
loadList = {
    '190310_I127_1_outfile.mat'
    '190415_I127_1_outfile.mat'
    '190418_I127_1_outfile.mat'
    '190420_I127_1_outfile.mat'
    '190429_I127_1_outfile.mat'
    '190501_I127_1_outfile.mat'
    '190506_I127_1_outfile.mat'
    '190513_I127_1_outfile.mat'
    '190515_I127_1_outfile.mat'
    '190523_HB40_1_outfile.mat'
    '190525_I127_1_outfile.mat'
    '190529_I127_1_outfile.mat'
    '190604_I127_1_outfile.mat'
    '190606_HB47_outfile.mat'
    '190826_I132_2_outfile.mat' %CamK2
    '191003_I135_1_outfile.mat' %camK2
    '191004_I135_1_outfile.mat' %CamK2z
    '191008_I135_1_outfile.mat' %CamK2
    '191014_I136_1_outfile.mat'
    '191015_W4_3_outfile.mat'
    '191016_I136_1_outfile.mat'
    '191017_I136_1_outfile.mat'
    '191121_I137_3_outfile.mat'
    '191122_I137_3_outfile.mat'
    '191126_I138_1_outfile.mat'
    '191206_I138_1_outfile.mat'
    '191212_I138_1_outfile.mat'
    '191213_W20_1_outfile.mat'
    '191217_mora_tre_outfile.mat'
    '191230_w20_1_outfile.mat'
    '200217_w14_1_outfile.mat'
    '200224_w14_1_outfile.mat'
    '200227_W14_1_outfile.mat'
    '200302_W14_1_outfile.mat'
    '200304_W14_1_outfile.mat'
    '200309_w21_1_outfile.mat'
    '200310_HB95_outfile.mat'
    '200311_i139_2_outfile.mat'
    '200312_i139_2_outfile.mat'
    '200316_i139_2_outfile.mat'
    '200728_i140_2_outfile.mat'
    '200723_i140_2_outfile.mat'
    '200729_I140_2_outfile.mat'
    '200807_w26_1_outfile.mat'
    '200810_i140_2_outfile.mat'
    '200902_w18_3_outfile.mat'
    '200901_w19_1_outfile.mat'
    '201102_w29_1_outfile.mat'
    '201103_w29_1_outfile.mat'
    ...'201106_w29_1_outfile.mat'% has 12 oris
    '201112_w29_3_outfile.mat'
    '201116_w29_3_outfile.mat'
    '201202_w29_3_outfile.mat'
    '210428_w32_2_outfile.mat' %used optimizer
    '210511_w32_2_outfile.mat' %many random holos for Sang Min not many reps might all be excluded
    '210830_w37_2_outfile.mat'
    '210902_I151_3_outfile.mat' %Many random holos at 30hz for SM Proj
    '210903_I151_3_outfile.mat' 
    '210826_i151_3_outfile.mat'
    '210831_w37_2_outfile.mat'
    '210901_w37_2_outfile.mat'
    ...'210910_i151_3_outfile.mat' %many random holos at 30Hz, for SM proj, (30c x 10ap)
    '210914_i151_3_outfile.mat'
    '210927_I154_2_outfile.mat'
    ...'210930_I154_2_outfile.mat' %has a weird not correctable error during respopnse window %exp2 doesn't suffer this
        '211021_W40_2_outfile.mat' %rate Exp, but several valid

    %manifold
   ...     '200722_W18_1_outfile.mat' %Spike Test Only (w Rates)
   ... '200724_W18_1_outfile.mat' %Spike Test Only
    ...'200723_I140_2_outfile.mat' %Spike Test Only (also 100spk)
    ...'201229_w29_4_outfile.mat' %spike Test Only

   ... '200730_I140_2_outfile.mat' %was this the example? %note there is a second mani attempt in this same recording, and a 'perceptual bias'
    ...'210105_W29_1_outfile.mat' %Note: 2 unexported manifold attempts
   ... '210111_W29_3_outfile.mat' %chrome2s Note: Mani is 75 Targets, Mani2 34 but all vis responsive
   ... '210127_W29_3_outfile.mat' %chrome2s
    %ai203 manifold
   ... '200727_W26_1_outfile.mat'
   ... '200728_W26_1_outfile.mat'
    ... '210426_I143_outfile.mat'
    ...'210517_I147_outfile.mat'
   ... '210518_I147_outfile.mat'
   ... '220127_HB120_outfile.mat'
   ... '220131_HB120_outfile.mat' %mani (daq epoch 10) is better used peg zero
   ... '220131_HB113_outfile.mat'
    
    %sepW1 load List
     '211108_I156_1_outfile.mat' %sepW1
    '211014_I156_1_outfile.mat'
    '211019_I156_1_outfile.mat'
    '211102_I158_1_outfile.mat'
    '211206_I162_1_outfile.mat' %very small dataset
    '211213_I162_1_outfile.mat' %Imaging looked weird, very suspect. 
    '220110_I162_2_outfile.mat' %sepW1 mostly 1 at a time 
    
    
    };


%% Deliberately Not Included
% '190307_I127_1_outfile.mat' %Too Different, stimparams don't exist right.
%    '190308_I127_1_outfile.mat' %Doesn't have a stimParam.roi and very
%    different format than usual
%    '190319_I127_1_outfile.mat' %Doesn't have a stimParam.roi and very
%    different format than usual
