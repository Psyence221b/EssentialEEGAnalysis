
function SPREP_func(subj_num,EEG,p,Purpose,alphaOn,Condition)
% Purpose 1 = sLORETA txt file output
% Purpose 2 = ERSP calculation
% Purpose 3 = ERP topofunction

if Purpose == 1
    for num_epoch = 1:size(EEG.data,3)
        Epoch = EEG.data(1:62,:,num_epoch)'; % ' is important!
        savename = sprintf('%s_CIT_%s_epoch %s.txt',subj_num,Condition,num2str(num_epoch));
        dlmwrite(savename,Epoch,'\t');
    end
    cd('..'); %go to Num \ sLORETA \ ERSP or ERP
    EEGave = mean(EEG.data,3)'; % be careful to ' <- (this is required to match loreta format)
    savename = sprintf('%s_CIT_%s_AVER.txt',subj_num,Condition);
    dlmwrite(savename,Epoch,'\t'); clear EEGave;
    
    
elseif Purpose == 2 % calculate ERSP 
if alphaOn == 1
    for i = 1:EEG(1).nbchan-4
    %calculate ERSP, ITC with bootstrap 0.05 confidence level   
    figure; [v_ERSP, v_ITC, v_powbase v_times v_frequencies] = pop_newtimef( EEG, 1, i, p.wavelet.time, p.wavelet.cycle , 'topovec', i, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', EEG.chanlocs(i).labels, 'baseline',p.wavelet.baseline,'alpha',p.wavelet.alpha,'freqs',p.wavelet.freqs,'nfreqs',p.wavelet.nfreqs, 'erspmax', p.wavelet.erspmax, 'itcmax',1, 'freqscale', 'log', 'plotphase', 'off', 'padratio', 1, 'winsize', 275);
    %save ERSP & ITC plotting
    savename = sprintf('[%s]ERSP_%s_%s.png',datestr(now,'yyyymmdd'),EEG.chanlocs(i).labels,Condition);
    saveas(gcf,savename); close figure 2;
    %save ERSP data
    savename = sprintf('[%s]ERSP_%s_%s',datestr(now,'yyyymmdd'),EEG.chanlocs(i).labels,Condition);
    cd('ERSP data'); save(savename, 'v_ERSP'); cd('..');
    %save ITC data
    savename = sprintf('[%s]ITC_%s_%s',datestr(now,'yyyymmdd'),EEG.chanlocs(i).labels,Condition);
    cd('ITC data'); save(savename, 'v_ITC'); cd('..');
    end
    cd('ERSP data'); save('Time','v_times'); save('Frequency','v_frequencies'); cd('..');
    cd('ITC data'); save('Time','v_times'); save('Frequency','v_frequencies'); cd('..');
else   
    for i = 1:EEG(1).nbchan-4
    %calculate ERSP, ITC    
    figure; [v_ERSP, v_ITC, v_powbase v_times v_frequencies] = pop_newtimef( EEG, 1, i, p.wavelet.time, p.wavelet.cycle , 'topovec', i, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', EEG.chanlocs(i).labels, 'baseline',p.wavelet.baseline,'freqs',p.wavelet.freqs,'nfreqs',p.wavelet.nfreqs, 'erspmax', p.wavelet.erspmax, 'itcmax',1, 'freqscale', 'log', 'plotphase', 'off', 'padratio', 1, 'winsize', 275);
    %save ERSP & ITC plotting
    savename = sprintf('[%s]ERSP_%s_%s.png',datestr(now,'yyyymmdd'),EEG.chanlocs(i).labels,Condition);
    saveas(gcf,savename); close figure 2;
    %save ERSP data
    savename = sprintf('[%s]ERSP_%s_%s',datestr(now,'yyyymmdd'),EEG.chanlocs(i).labels,Condition);
    cd('ERSP data'); save(savename, 'v_ERSP'); cd('..');
    %save ITC data
    savename = sprintf('[%s]ITC_%s_%s',datestr(now,'yyyymmdd'),EEG.chanlocs(i).labels,Condition);
    cd('ITC data'); save(savename, 'v_ITC'); cd('..');
    end
    %save times, frequencies variables
    cd('ERSP data'); save('Time','v_times'); save('Frequency','v_frequencies'); cd('..');
    cd('ITC data'); save('Time','v_times'); save('Frequency','v_frequencies'); cd('..');
end

elseif Purpose == 3 % ERP topofunction
    EEG = pop_select( EEG,'nochannel',{'UVEO' 'LVEO' 'RHEO' 'LHEO'});
    EEG = eeg_checkset( EEG ); eeglab redraw;
    pop_topoplot(EEG,1, [-200:100:1000] ,'Target',[1 13] ,0,'electrodes','off','maplimits',[-12 12] );
    savename = sprintf('[%s]Topo ERP_%s.png',datestr(now,'yyyymmdd'),Condition);
    pic = getframe(gcf);
    imwrite(pic.cdata,savename); close figure 2; 
end



end


