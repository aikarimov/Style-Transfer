FName = {'VanGogh_Chestnut','VanGogh_Wheat','VanGogh_Irises','Seurat_Chahut','Seurat_Bridge','Levitan_Oak','Levitan_Evening','Lenna'};
nfiles = length(FName);
ny = 2; nx = 4; 
for n=1:8
    Im = imread([FName{1,n},'.jpg']);
    figure(11); %нарисуем гистограмму длин
    subplot(ny,nx,n);
    imshow(Im);
    title(['(','a' + n - 1,')']);
end