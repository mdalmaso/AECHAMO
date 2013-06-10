function subplot_dmps(obj,sub,varargin)
% SUBPLOT_DMPS plots the distribution.

% (c) Miikka Dal Maso 2013
%
% Version history:
% 2013-05-24    0.1.0
% 2013-06-10    0.1.1 Takes now argument 'smoothed' to plot the smoothed
%                     distribution instead of the original.

if(nargin > 2)
    if(strcmp(varargin{1},'smoothed'))
        v = obj.output_data.distr_smoothed;
    elseif(nargin > 3)
        error('Too many arguments.');
    else
        error('Invalid argument: ''%s''.',varargin{1});
    end
else
    v = obj.output_data.distr;
end

%% first one
[m n]=size(v);

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

Zdata=log10(abs(v(2:m,3:n))+1e-21);

pcolor(v(2:m,1),v(1,3:n),Zdata');
% shading interp

 colormap(jet(250));

%rr = 0:0.15:1;
%rr = flipud(rr');

%colm = [rr rr rr];
%colormap(colm);

caxis([1 4]);

shading flat
% rr = 0:0.15:1;
% rr = flipud(rr');

% colm = [rr rr rr];




% colormap(colm);

% caxis([1 5]);


set(gca,'yscale','log')
% axis([fix(V(:,1)) fix(V(:,1))+1 1e-9 1e-6]);
% set(gca,'xtick',xti,'xticklabel',xtila,'fontsize',14)
grid


% gcc=colorbar('horiz');
% set(gcc,'xlim',[1 5],'xtick',[1 2 3 4 5],'fontsize',14)
% set(gcc,'xticklabel',[10,100,1000,10000 100000]')
% V=get(gcc,'position');
% set(gcc,'position',[V(1) V(2)-0.03 V(3) V(4)])



%% first one



