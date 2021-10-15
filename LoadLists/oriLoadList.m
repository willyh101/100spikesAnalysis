% add experiments for ori ensembles here...
loadListName = 'oriLoadList'; 

loadList = {
    '190418_I127_1_outfile.mat'
    '190506_I127_1_outfile.mat'
    '190513_I127_1_outfile.mat'
    '190515_I127_1_outfile.mat'
    '190523_HB40_1_outfile.mat'
    '190525_I127_1_outfile.mat'
    '190529_I127_1_outfile.mat'
    '190604_I127_1_outfile.mat'
    ...'190606_HB47_outfile.mat' %ai203 (exclude)
    '191017_I136_1_outfile.mat'
    '191126_I138_1_outfile.mat'
    '191206_I138_1_outfile.mat'
    '191212_I138_1_outfile.mat'
    ...'191213_W20_1_outfile.mat' %ai203 (exclude)
    '191217_mora_tre_outfile.mat'
    ...'191230_w20_1_outfile.mat' %ai203 (exclude)
    '200302_W14_1_outfile.mat' 
    '200304_W14_1_outfile.mat' 
    '200309_w21_1_outfile.mat'
    ...'200310_HB95_outfile.mat'  %ai203 (exclude)
    '200311_i139_2_outfile.mat' %will excluded
    '200728_i140_2_outfile.mat' 
    '200723_i140_2_outfile.mat'
    '200729_I140_2_outfile.mat'
    ...'200807_w26_1_outfile.mat'  %ai203 (exclude)
    '200810_i140_2_outfile.mat'
    '200902_w18_3_outfile.mat'
    '200901_w19_1_outfile.mat'
    '201102_w29_1_outfile.mat'
    '201103_w29_1_outfile.mat'
%     '201106_w29_1_outfile.mat'% 12 oris
    '201112_w29_3_outfile.mat'
    '201116_w29_3_outfile.mat'
    '201202_w29_3_outfile.mat'
    '210428_w32_2_outfile.mat' % first of optimizer
    '210830_w37_2_outfile.mat' %neo-IV Tre2s
    '210902_I151_3_outfile.mat' %Many random holos at 30hz for SM Proj (strangely many failure..?)
    '210903_I151_3_outfile.mat' 
    '210826_i151_3_outfile.mat'
    '210831_w37_2_outfile.mat' %neo-IV Tre2s
    '210901_w37_2_outfile.mat' %neo-IV Tre2s
    ...'210910_i151_3_outfile.mat' %many random holos at 30Hz, for SM proj, (30c x 10ap)
    '210914_i151_3_outfile.mat'
    '210927_I154_2_outfile.mat'
    ...'210930_I154_2_outfile.mat' %has a weird not correctable error during respopnse window
};