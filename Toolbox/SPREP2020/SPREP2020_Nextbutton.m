function Nextbutton( Name, String )
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name',Name,'NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String',String,...
          'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check

end

