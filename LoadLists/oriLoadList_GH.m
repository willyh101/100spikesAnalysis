function [loadList] = oriLoadList_GH(listType)

% Smaller load list to debug the code
if strcmp(listType,'short')
    loadList = {'210428_w32_2_outfile.mat'
    '210903_I151_3_outfile.mat'
    '210826_i151_3_outfile.mat'
    '210914_i151_3_outfile.mat'
    '210927_I154_2_outfile.mat'
    '211019_I156_1_outfile.mat'};

elseif strcmp(listType,'long')
    loadList = {'190418_I127_1_outfile.mat'
    '190506_I127_1_outfile.mat'
    '190513_I127_1_outfile.mat'
    '190515_I127_1_outfile.mat'
    '190525_I127_1_outfile.mat'
    '190529_I127_1_outfile.mat'
    '190604_I127_1_outfile.mat'
%     '191017_I136_1_outfile.mat' % something off with the # of time frames in tracesVis
%     '191126_I138_1_outfile.mat' %rate exp, but several valid; vis stim and holo stim off by 1 total frame
    '191206_I138_1_outfile.mat'
    '191217_mora_tre_outfile.mat'
    '200302_W14_1_outfile.mat'
    '200304_W14_1_outfile.mat'
    '200309_w21_1_outfile.mat'
    '200311_i139_2_outfile.mat'
%     '200728_i140_2_outfile.mat' % vis stim and holo stim off by 1 total frame (31 vs. 30)
    '200723_i140_2_outfile.mat'
    '200729_I140_2_outfile.mat'
    '200810_i140_2_outfile.mat'
    '200902_w18_3_outfile.mat'
    '200901_w19_1_outfile.mat'
    '201103_w29_1_outfile.mat'
%    '201112_w29_3_outfile.mat' % vis stim and holo stim off by 1 total frame
    '201116_w29_3_outfile.mat'
%    '201202_w29_3_outfile.mat' % vis stim and holo stim off by 1 total frame
    '210428_w32_2_outfile.mat' % first of optimizer
    '210902_I151_3_outfile.mat' %Many random holos at 30hz for SM Proj (strangely many failure..?)
    '210903_I151_3_outfile.mat'
    '210826_i151_3_outfile.mat'
    '210914_i151_3_outfile.mat'
    '210927_I154_2_outfile.mat'
    '211021_W40_2_outfile.mat' %rate Exp, but several valid
    '211019_I156_1_outfile.mat' %sepW1
    '211102_I158_1_outfile.mat' %sepW1
    '211108_I156_1_outfile.mat' %sepW1  one at a time
    '211019_I154_1_outfile.mat' %one at a time but some eligible
    };
else
    error('Unknown loadList. Options are ''short'' or ''long'' ')
end

