function subplot_dmps(obj,sub,varargin)
% SUBPLOT_DMPS plots the distribution.

% (c) Miikka Dal Maso 2013
%
% Version history:
% 2013-05-24    0.1.0
% 2013-06-10    0.1.1 Takes now argument 'smoothed' to plot the smoothed
%                     distribution instead of the original.

% if(nargin > 2)
%     if(strcmp(varargin{1},'smoothed'))
%         v = obj.output_data.distr_smoothed;
%     elseif(nargin > 3)
%         error('Too many arguments.');
%     else
%         error('Invalid argument: ''%s''.',varargin{1});
%     end
% else
%     v = obj.output_data.distr;
% end

plot_original = 0;

if(nargin > 2)
    for i=1:nargin-2
        switch(varargin{i})
            case('smoothed')
                v = obj.output_data.distr_smoothed;
                Zdata=log10(abs(v(2:m,3:n))+1e-21);
            case('original')
                v = obj.output_data.distr_original;
                Zdata=log10(abs(v(2:end,3:(end-2)/2+2))+(1e-21));
                for i=1:length(Zdata(1,:))
                    time(:,i)=v(2:end,1);
                end
                plot_original = 1;
            otherwise
                error('Invalid argument: ''%s''.',varargin{i});
        end
    end
else
    v = obj.output_data.distr;  
    Zdata=log10(abs(v(2:end,3:end))+1e-21);
end


%% first one
% [m n]=size(v);

%% calculating <15 nm particles

% V = Dlog_to_N(v);



% fprintf('\n%2i',0)
% for i = 2:m,
% fprintf('\b\b\b%2i',i)
% pause(0.1)
% end;

figure(gcf)

%set(gcf,'papertype','a4letter')

subplot(sub)

% set time from seconds to hours
v(2:end,1) = v(2:end,1)/3600;

%set time from hours to days
v(2:end,1) = v(2:end,1)/24;

semilogy(v(2:end,1),v(2:end,2),'or');
V=axis;
ve=fix(V(:,1));
%xti=[ve:0.125:ve+1];

for ii=1:2:9,
   lab=num2str(ve+(ii-1)*0.125,5);
   xtila(ii,1:length(lab))=lab;
end;

%set(gca,'xtick',xti,'xticklabel',xtila,'fontsize',14)
grid


subplot(sub)


if(plot_original == 1)
    pcolor(time,v(2:end,(end-2)/2+3:end),Zdata);
else
    pcolor(v(2:end,1),v(1,3:end),Zdata');
end
% shading interp

 colormap(jet(250));

%rr = 0:0.15:1;
%rr = flipud(rr');

%colm = [rr rr rr];
%colormap(colm);

%caxis([1 4.5]);
%colorbar

shading flat
% rr = 0:0.15:1;
% rr = flipud(rr');

% colm = [rr rr rr];




% colormap(colm);

 caxis([1 5]);


set(gca,'yscale','log')
% axis([fix(V(:,1)) fix(V(:,1))+1 1e-9 1e-6]);
% set(gca,'xtick',xti,'xticklabel',xtila,'fontsize',14)
grid

%% colorbar fix
gcc=colorbar('horiz');
set(gcc,'xlim',[1 5],'xtick',[1 2 3 4 5],'fontsize',10)
set(gcc,'xticklabel',[10,100,1000,10000 100000]')
x = xlabel(gcc, 'dN/dlogDp (cm^{-3})');
set(x, 'fontsize',12)
%set(x,'Fontangle','italic')

% V=get(gcc,'position');
% set(gcc,'position',[V(1) V(2)+0.83 V(3) V(4)])

%% vertical colorbar
% gcc=colorbar('vertic');
% set(gcc,'ylim',[1 5],'ytick',[1 2 3 4 5],'fontsize',14)
% set(gcc,'yticklabel',[10,100,1000,10000 100000]')
% % V=get(gcc,'position');
% % set(gcc,'position',[V(1) V(2)-0.01 V(3) V(4)])

%% x-axis fix
xhandle = xlabel('time (d)');
set(xhandle,'Fontsize',12)
%set(xhandle,'Fontangle','italic')
%set(xhandle,'Fontname','Computermodern')

%% y-axis fix
yhandle = ylabel('D_{p} (nm)','rotation',90);
set(yhandle,'Fontsize',12)
%set(yhandle,'Fontangle','italic')
% set dp from m to nm
set(gca,'YTickLabel',[1,10,100, 1000])
%set(yhandle,'Fontname','Computermodern')

%% first one

%% Testiä
% v=[dist(2:end,1), dist_original(:,:)];
% Zdata=log10(abs(v(1:end,2:end))+1e-21);
% pcolor(v(1:end,1),v(1:end,41:end),Zdata');
% Error using pcolor (line 57)
% Matrix dimensions must agree.
%  
% Zdata=log10(abs(v(1:end,2:41))+1e-21);
% pcolor(v(1:end,1),v(1:end,41:end),Zdata');
% Error using pcolor (line 57)
% Matrix dimensions must agree.
%  
% length(v(1,41:end))
% 
% ans =
% 
%     41
% 
% pcolor(v(1:end,1),v(1:end,42:end),Zdata');
% Error using pcolor (line 57)
% Matrix dimensions must agree.
%  
% length(v(1,42:end))
% 
% ans =
% 
%     40
% 
% Zdata=log10(abs(v(1:end,2:41))+(1e-21));
% 1e-21
% 
% ans =
% 
%    1.0000e-21
% 
% size(Zdata')
% 
% ans =
% 
%     40   301
% 
% size(v(1:end,2:41))
% 
% ans =
% 
%    301    40
% 
% pcolor(v(1:end,1),v(1:end,42:end),Zdata);
% Error using pcolor (line 57)
% Matrix dimensions must agree.
%  
% timemat=v(1:end,1).*ones(size(Zdata));
% Error using  .* 
% Matrix dimensions must agree.
%  
% for i=1:length(Zdata(1,:))
% timemat(:,i)=v(1:end,1);
% end



