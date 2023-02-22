clear all;useCase='Kallerup';dx=0.25;forward.type='ray';is_slowness=0;caseTomo_setup
clear all;fmat='caseTomo_Kallerup_dx25_Fray-none_ME0_slo0.mat';N=500000;di_use=1;caseTomo_metropolis
%prior_reals=1./prior_reals;
%%


%%
nx=length(prior{1}.x);
ny=length(prior{1}.y);
for i=1:size(reals_all,1);
    m{1}=reshape(reals_all(i,:),ny,nx);
    d=sippi_forward(m,forward,prior,data);
    d_post(:,i)=d{1};
    dd(:,i)=data{1}.d_obs-d{1};
end

d_std=sqrt(diag(data{1}.Ct))
dd_rel = std(dd')'./d_std;

figure(1);
subplot(3,1,1)
plot(d_std)
hold on
plot(std(dd'))
hold off
ylabel('Std')
legend('Noise','Post')

subplot(3,1,2)
plot(dd_rel)
ylabel('RelStd')

subplot(3,1,3)
plot(dd_rel)
plot(d_post,'r-');
hold on
plot(data{1}.d_obs,'k-','LineWidth',2);
hold off
grid on
ylabel('Traveltime (ns)')
xlabel('Data #')



%%
figure(9);clf
plot(forward.sources(:,1),forward.sources(:,2),'k.')
hold on
plot(forward.receivers(:,1),forward.receivers(:,2),'k.')
cmap=flipud(hot);
cax=[min(dd_rel) max(dd_rel)];
colormap(gca,cmap);caxis(cax);
cb=colorbar;
set(get(cb,'Ylabel'),'String','REL ERR')
vapp=dd_rel
for i=1:size(forward.sources,1)
    %icol=ceil(interp1([cax],[1 size(cmap,1)],vapp(i),'nearest'))
    icol(i)=ceil(interp1([cax],[1 size(cmap,1)],vapp(i),'linear','extrap'));
    if vapp(i)<cax(1);icol(i)=1;end
    if vapp(i)>cax(2);icol(i)=size(cmap,1);end
    lw=2*(vapp-0.1)./(0.1);
    lw=5*(vapp-0.4)./(0.8)
    lw(lw<0.01)=0.01;

    plot([forward.sources(i,1),forward.receivers(i,1)],[forward.sources(i,2),forward.receivers(i,2)],'-','LineWidth',lw(i),'MarkerSize',1,'Color',cmap(icol(i),:))
    %plot([D.S(i,1),D.R(i,1)],[D.S(i,2),D.R(i,2)],'-','LineWidth',.1)
end
grid on
hold off
axis image
set(gca,'ydir','reverse')
print_mul(sprintf('%s_%s',txt,'d_Err'))




