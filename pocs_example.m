clear
close all
N = 256;
iter = 20;
pf = 9/16;
% Set up some pseudo coilmaps:
cmap(1,:,:) = repmat(        ( 3./(1:N)  .^0.10) .* exp(1i*(1:N)  *pi/(N/3)) , N, 1 );
cmap(2,:,:) = repmat(        ( 1./(1:N).'.^0.20) .* exp(1i*(1:N).'*pi/(N/4)) , 1, N );
cmap(3,:,:) = repmat( fliplr(( 2./(1:N)  .^0.04) .* exp(1i*(1:N)  *pi/(N/8))), N, 1 );
cmap(4,:,:) = repmat( flipud(( 2./(1:N).'.^0.15) .* exp(1i*(1:N).'*pi/(N/10))),1, N );
cmap = double(cmap);

Nc = size(cmap,1);
p = reshape( double(phantom(N)), 1, N, N );
subs = { 1, 164:181, 187:193 };
p(subs{:}) = 0;
subs = { 1, 166:179, 189:191 };
p(subs{:}) = 1 * exp(1i * pi/2);
p = bsxfun( @times, cmap, p );
kspRef = fftshift(fftshift(  fft(fft(  fftshift(fftshift( p, 2),3), [],2),[],3), 2),3);
kspRef = kspRef + 10 * complex( randn(Nc,N,N), randn(Nc,N,N) );    % add noise
kspPF  = kspRef;
kspPF(:,:,1:floor((1-pf)*N)) = 0;
sharpness = zeros(iter, 1);
imReco_iter = zeros(iter, N, N);

for i = 1:iter
    imReco = pocs(kspPF, i, 0^logical(iter-i));
    imReco   = squeeze(sqrt(sum( abs(imReco).^2, 1)));
    imReco_iter(i, :, :) = squeeze(imReco);
    sharpness(i) = calculateSharpness(imReco);
end
imRef    = ifftshift(ifftshift(  ifft(ifft(  ifftshift(ifftshift( kspRef, 2),3), [],2),[],3), 2),3);
imRef    = squeeze(sqrt(sum( abs(imRef).^2, 1)));
imDirect = ifftshift(ifftshift(  ifft(ifft(  ifftshift(ifftshift(  kspPF, 2),3), [],2),[],3), 2),3);
imDirect = squeeze(sqrt(sum( abs(imDirect).^2, 1)));

figure;
plot(1:iter, sharpness, 'o-');
xlim([1 iter])
xlabel('Iteration');
ylabel('Sharpness');
title('Sharpness Improvement');
set(gca,"FontName","Times")

figure,
imagesc(abs([imRef  imReco  imDirect]),[0 2*mean(abs(imRef(:)))])
%imshow
colorbar
axis image
title('reference (full dataset)   |   POCS reco   |    zero-padded reco')

figure,
imshow(abs([imReco-imRef  imDirect-imRef]))
axis image
title('POCS reco   |    zero-padded reco')

disp(["reference (full dataset)","POCS reco","zero-padded reco";...
      calculateSharpness(imRef),calculateSharpness(imReco),calculateSharpness(imDirect)])

%% Extracting differences in images by filtering
block_size = 3;  % Size of the block (3x3 for 9 pixels)
X = 1:block_size:size(imRef, 1);
Y = 1:block_size:size(imRef, 2);
if block_size>1
X(X>=256)=[];
Y(Y>=256)=[];
end
diff_values_mean = zeros([length(X),length(Y),2]);
diff_values_min = zeros([length(X),length(Y),2]);
diff_values_max = zeros([length(X),length(Y),2]);

ii = 0;
jj = 0;
for i = X
    ii = ii + 1;
    jj = 0;
    for j = Y
        jj = jj + 1;
        block_Reco = double(imReco(i:i+block_size-1, j:j+block_size-1));
        block_Direct = double(imDirect(i:i+block_size-1, j:j+block_size-1));
        block_Ref = double(imRef(i:i+block_size-1, j:j+block_size-1));
        diff_block_Reco = abs(block_Reco - block_Ref);
        diff_block_Direct = abs(block_Direct - block_Ref);
        diff_value_mean = [mean(diff_block_Reco(:));mean(diff_block_Direct(:))];
        diff_value_min = [min(diff_block_Reco(:));min(diff_block_Direct(:))];
        diff_value_max = [max(diff_block_Reco(:));max(diff_block_Direct(:))];
        diff_values_mean(ii,jj,:) = diff_value_mean;
        diff_values_min(ii,jj,:) = diff_value_min;
        diff_values_max(ii,jj,:) = diff_value_max;
    end
end
figure,
imshow([diff_values_mean(:,:,1),diff_values_mean(:,:,2);...
        diff_values_min(:,:,1),diff_values_min(:,:,2);...
        diff_values_max(:,:,1),diff_values_max(:,:,2)])
axis image
title('POCS reco   |    zero-padded reco')
ylabel('max   |   min   |   mean')

figure,
yyaxis left
semilogy(X,abs(mean(diff_values_mean(:,:,1))),DisplayName='POCS',LineWidth=1.5);hold on
plot(X,abs(mean(diff_values_mean(:,:,2))),DisplayName='Direct',LineWidth=1);hold off
xlim([min(X) max(X)])
ylim([min(min(abs(mean(diff_values_mean)))) 1.1*max(max(abs(mean(diff_values_mean))))])
xlabel('Pixels')
ylabel('Average of the difference in the x direction')
yyaxis right
semilogy(Y,abs(mean(diff_values_mean(:,:,1),2)),DisplayName='POCS',LineWidth=1.5);hold on
plot(Y,abs(mean(diff_values_mean(:,:,2),2)),DisplayName='Direct',LineWidth=1);xlim([min(Y) max(Y)])
ylim([min(min(abs(mean(diff_values_mean,2)))) 1.1*max(max(abs(mean(diff_values_mean,2))))])
ylabel('Average of the difference in the y direction')
legend
title(['Mean error of filtering by ',num2str(block_size),'\times',num2str(block_size),'pixels'])
set(gca,FontName="Times")

figure,
yyaxis left
semilogy(X,abs(mean(diff_values_min(:,:,1))),DisplayName='POCS',LineWidth=1.5);hold on
plot(X,abs(mean(diff_values_min(:,:,2))),DisplayName='Direct',LineWidth=1);hold off
xlim([min(X) max(X)])
ylim([min(min(abs(mean(diff_values_min)))) 1.1*max(max(abs(mean(diff_values_min))))])
xlabel('Pixels')
ylabel('Average of the difference in the x direction')
yyaxis right
semilogy(Y,abs(mean(diff_values_min(:,:,1),2)),DisplayName='POCS',LineWidth=1.5);hold on
plot(Y,abs(mean(diff_values_min(:,:,2),2)),DisplayName='Direct',LineWidth=1);xlim([min(Y) max(Y)])
ylim([min(min(abs(mean(diff_values_min,2)))) 1.1*max(max(abs(mean(diff_values_min,2))))])
ylabel('Average of the difference in the y direction')
legend
title(['Minimum error of filtering by ',num2str(block_size),'\times',num2str(block_size),'pixels'])
set(gca,FontName="Times")

figure;
yyaxis left
semilogy(X,abs(mean(diff_values_max(:,:,1))),DisplayName='POCS',LineWidth=1.5);hold on
plot(X,abs(mean(diff_values_max(:,:,2))),DisplayName='Direct',LineWidth=1);hold off
xlim([min(X) max(X)])
ylim([min(min(abs(mean(diff_values_max)))) 1.1*max(max(abs(mean(diff_values_max))))])
xlabel('Pixels')
ylabel('Average of the difference in the x direction')
yyaxis right
semilogy(Y,abs(mean(diff_values_max(:,:,1),2)),DisplayName='POCS',LineWidth=1.5);hold on
plot(Y,abs(mean(diff_values_max(:,:,2),2)),DisplayName='Direct',LineWidth=1);xlim([min(Y) max(Y)])
ylim([min(min(abs(mean(diff_values_max,2)))) 1.1*max(max(abs(mean(diff_values_max,2))))])
ylabel('Average of the difference in the y direction')
legend
title(['Maximum error of filtering by ',num2str(block_size),'\times',num2str(block_size),'pixels'])
set(gca,FontName="Times")

function sharpness = calculateSharpness(image)
    [gradX, gradY] = gradient(double(image));
    sharpness = mean(abs(gradX(:))) + mean(abs(gradY(:)));
end