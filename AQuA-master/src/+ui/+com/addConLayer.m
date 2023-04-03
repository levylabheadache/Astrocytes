function addConLayer(f,bLayer)
% movie
bL1 = uix.VBox('Parent',bLayer,'Spacing',1,'Padding',3);
uicontrol(bL1,'Style','text','String','--- Movie brightness/contrast ---');
% uix.Empty('Parent',bL1);
uicontrol(bL1,'Style','text','String','Min','HorizontalAlignment','left');
uicontrol(bL1,'Style','slider','Tag','sldMin','Callback',{@ui.over.adjMov,f});
uix.Empty('Parent',bL1);
uicontrol(bL1,'Style','text','String','Max','HorizontalAlignment','left');
uicontrol(bL1,'Style','slider','Tag','sldMax','Callback',{@ui.over.adjMov,f});
uix.Empty('Parent',bL1);

% single view
pBrightness = uix.CardPanel('Parent',bL1,'Tag','pBrightness');  % two sliders
pB = uix.VBox('Parent',pBrightness);
uicontrol(pB,'Style','text','String','Brightness','HorizontalAlignment','left');
uicontrol(pB,'Style','slider','Tag','sldBri','Callback',{@ui.over.adjMov,f});
pB.Heights = [15,15];

% side by side view
pBS = uix.Grid('Parent',pBrightness);
uicontrol(pBS,'Style','text','String','Left Brightness','HorizontalAlignment','left');
uicontrol(pBS,'Style','slider','Tag','sldBriL','Callback',{@ui.over.adjMov,f});
uix.Empty('Parent',pBS);
uix.Empty('Parent',pBS);
uicontrol(pBS,'Style','text','String','Right Brightness','HorizontalAlignment','left');
uicontrol(pBS,'Style','slider','Tag','sldBriR','Callback',{@ui.over.adjMov,f});
pBS.Widths = [-1,2,-1];   pBS.Heights = [15,15];

pBrightness.Selection = 1;
bL1.Heights = [18,15,15,3,15,15,3,30];

uix.Empty('Parent',bLayer);

% overlays
bL2 = uix.VBox('Parent',bLayer,'Spacing',1,'Padding',3);
uicontrol(bL2,'Style','text','String','--- Feature overlay ---');
% uix.Empty('Parent',bL2);
uicontrol(bL2,'Style','text','String','Type','HorizontalAlignment','left');
uicontrol(bL2,'Style','popupmenu','Tag','overlayDat','String',{'None'},'Callback',{@ui.over.chgOv,f,1});
uix.Empty('Parent',bL2);
x0 = [18,15,20,3];

uicontrol(bL2,'Style','text','String','Feature','HorizontalAlignment','left');
uicontrol(bL2,'Style','popupmenu','Tag','overlayFeature','String',{'Index'},'Enable','off');
uicontrol(bL2,'Style','text','String','Color','HorizontalAlignment','left');
uicontrol(bL2,'Style','popupmenu','Tag','overlayColor','String',{'Random','GreenRed'},'Enable','off');
x1a = [15,20,15,20];
uicontrol(bL2,'Style','text','String','Transform','HorizontalAlignment','left');
uicontrol(bL2,'Style','popupmenu','Tag','overlayTrans','String',{'None','Square root','Log10'},'Enable','off');
uicontrol(bL2,'Style','text','String','Divide','HorizontalAlignment','left');
uicontrol(bL2,'Style','popupmenu','Tag','overlayScale','String',{'None','Size','SqrtSize'},'Enable','off');
x1b = [15,20,15,20];
uicontrol(bL2,'Style','text','String','Propagation direction','HorizontalAlignment','left');
uicontrol(bL2,'Style','popupmenu','Tag','overlayPropDi','String',...
    {'Anterior','Posterior','Lateral Left','Lateral Right'},'Enable','off');
uicontrol(bL2,'Style','text','String','Landmark ID','HorizontalAlignment','left');
uicontrol(bL2,'Style','edit','Tag','overlayLmk','String','1','Enable','off','HorizontalAlignment','left');
x1c = [15,20,15,20];
uix.Empty('Parent',bL2);
x1 = [x1a,x1b,x1c,5];

bDrawBt = uix.HButtonBox('Parent',bL2,'Spacing',5,'ButtonSize',[120,20]);
uicontrol(bDrawBt,'String','Update overlay','Tag','updtFeature','Callback',{@ui.over.chgOv,f,2});
uicontrol(bDrawBt,'String','Read features','Tag','addUserFeature','Callback',{@ui.over.chgOv,f,0});
uicontrol(bDrawBt,'String','Read colormaps','Callback',{@ui.over.getColMap,f});
uix.Empty('Parent',bL2);
x2 = [20,12];

uicontrol(bL2,'Style','text','String','Min','HorizontalAlignment','left','Tag','txtMinOv');
uicontrol(bL2,'Style','slider','Tag','sldMinOv','Callback',{@ui.over.adjMov,f,1},'Enable','off');
uicontrol(bL2,'Style','text','String','Max','HorizontalAlignment','left','Tag','txtMaxOv');
uicontrol(bL2,'Style','slider','Tag','sldMaxOv','Callback',{@ui.over.adjMov,f,1},'Enable','off');
uicontrol(bL2,'Style','text','String','Brightness','HorizontalAlignment','left');
uicontrol(bL2,'Style','slider','Tag','sldBriOv','Callback',{@ui.over.adjMov,f,1});
x3 = [15,15,15,15,15,15];

bL2.Heights = [x0,x1,x2,x3];
bLayer.Heights = [130,-1,440];
end





