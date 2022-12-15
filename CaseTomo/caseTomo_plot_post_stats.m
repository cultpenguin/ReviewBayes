function caseTomo_plot_post_stats(post_reals,prior_reals,prior,txt_out)

if nargin<4, txt_out = 'stat';end
if nargin<3, prior{1}.null=[];end
if ~isfield(prior{1},'cax');prior{1}.cax=[0.1150    0.1750];end
if ~isfield(prior{1},'cax_std');prior{1}.cax_std=[0 0.02];end
if ~isfield(prior{1},'cmap');prior{1}.cmap=jet;end

x=prior{1}.x;
y=prior{1}.y;
nx=length(prior{1}.x);
ny=length(prior{1}.y);
dx=prior{1}.x(2)-prior{1}.x(1);

Nr=5;
 

% plot posterior stats
figure(11);
for i=1:Nr
    subplot(1,Nr,i)
    imagesc(x,y,reshape(prior_reals(i,:),ny,nx))
    axis image
    caxis(prior{1}.cax);colormap(prior{1}.cmap)
    title('\rho(m)\rightarrowm^*')
end
colorbar_shift;
print_mul(sprintf('%s_prior_sample',txt_out))

figure(12);
i_plot = ceil(linspace(1,size(post_reals,1),Nr));
for i=1:Nr
    subplot(1,Nr,i)
    m=sippi_prior(prior);
    imagesc(x,y,reshape(post_reals(i_plot(i),:),ny,nx))
    axis image;caxis(prior{1}.cax);colormap(prior{1}.cmap)
    title('\sigma(m)\rightarrowm^*')
end
colorbar_shift;
print_mul(sprintf('%s_post_sample',txt_out))


%%
figure(13);clf
subplot(1,3,1)
imagesc(prior{1}.x,prior{1}.y,reshape(mean(post_reals),ny,nx))
axis image;caxis(prior{1}.cax);colormap(prior{1}.cmap)
colorbar
title('\sigma(m) - mean')
subplot(1,3,2)
imagesc(prior{1}.x,prior{1}.y,reshape(std(post_reals),ny,nx))
axis image;colormap(prior{1}.cmap)
caxis(prior{1}.cax_std)
title('\sigma(m) - std')
colorbar
print_mul(sprintf('%s_post_mean_std',txt_out))
    

%%
vmax=0.13;
P=mean(post_reals<vmax);
subplot(1,3,3)
imagesc(prior{1}.x,prior{1}.y,reshape(P,ny,nx))
axis image;colormap(gca,flipud(hot))
title(sprintf('P_{ \\sigma }(m<%3.2f)',vmax))
caxis([0 1])
colorbar
print_mul(sprintf('%s_post_mean_std_P',txt_out))

%% Area where V<vmax
Apixel=(post_reals<vmax).*(dx*dx);
A=sum(Apixel,2);
Amedian=median(A);
figure(14);clf
wb=2;
Abin=[4.5:wb:40.5];
Ac=(Abin(2:end)+Abin(1:end-1))/2;
h=histogram(A,Abin);
Apdf=h.Values;
Apdf=Apdf/(wb*sum(Apdf));
bar(Ac,Apdf)
xlabel(sprintf('Area where m<%3.2f',vmax))
ylabel('Probability')
hold on
yl=ylim;xl=xlim;
plot([1 1].*Amedian,ylim,'k:','LineWidth',3)
text(Amedian+0.02*diff(xl),yl(1)+0.95*diff(yl),sprintf('A_{median}=%3.1f',Amedian),'FontSize',10)
hold off
grid on
ppp(8,8,9,2,2)
print_mul(sprintf('%s_area',txt_out))


