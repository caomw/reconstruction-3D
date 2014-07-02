clear all;
rowCorners= zeros(19,1000) ;
corners=zeros(19,768,1024) ;
colCorners= zeros(19,1000) ;
projMatrix = zeros(19,4,4) ;

I= zeros(19,768,1024) ;
Ib = zeros(19,768,1024) ;

for i=1:1:19
    
    I(i,:,:) = imread(strcat('40/lybe_Sil_',int2str(i),'_40.pgm'));
    
    corners(i,:,:) = imregionalmax(cornermetric(reshape(I(i,:,:),768,1024)));
    
    [row2,col2] = find(reshape(corners(i,:,:),768,1024));
    rowCorners(i,1:size(row2,1)) = row2' ;
    colCorners(i,1:size(col2,1)) = col2' ;
end


for i= 1:1:19
    Ib(i,:,:) = imregionalmax( reshape(I(i,:,:),768,1024) ) ;
end

allMatrixes = load('calibfile');

for i=1:1:19
    projMatrix(i,1:3,1:4) = allMatrixes((3*(i-1)+1):(i*3),1:4) ;
    projMatrix(i,4,4) = 1 ;
end

cameraCenters = zeros(4,19);

for i = 1:19
    CC = reshape(projMatrix(i,1:3,:), 3, 4) ;
    [U, D, V] = svd(CC);
    cameraCenters(:, i) = V(:,end);
end

cameraCenters = cameraCenters ./ repmat(cameraCenters(4,:), 4,1);


pImage = ones([4 1]) ;
pImage2 = ones([4 1]) ;

fImage = ones(1,3) ;
fImage2 = ones(1,3) ;

pTemp = zeros(4,1) ;
pTemp2 = zeros(4,1) ;
pReal = zeros(4,1) ;
pReal2 = zeros(4,1) ;

Xs = zeros(2,600) ;
Ys = zeros(2,600) ;
Zs = zeros(2,600) ;

counter1 = 0 ;
coun = 0 ;
Xw = zeros(4,0);
yys = zeros( size(rowCorners,2), 1024 ) ;
yys2 = zeros( size(rowCorners,2), 1024 ) ;
intStarts = zeros(2,0) ;
intEnds = zeros(2,0) ;
triangX1 = zeros(2,0) ;
triangX2 = zeros(2,0) ;

lastXw = zeros(3,0) ;
for j= 1:1:19
    cameraNum = j ;
    FoundPointsNum = size(lastXw,2) ;
    for i=1:1:size(rowCorners,2)
        lastIntervals = zeros ( 6, 0 ) ;
        camCount = 0 ;
        for k = [(mod(j-5,19)+1) (mod(j-3,19)+1) (mod(j+1,19)+1) (mod(j+3,19)+1)]
            camCount = camCount + 1 ;
            if( camCount < 3)
                coef=-1 ;
            else
                coef=1 ;
            end
            
            pr1 = reshape(projMatrix(j,1:3,:),3,4) ;
            pr2 = reshape(projMatrix(k,1:3,:),3,4) ;
            
            if( rowCorners(j,i) ~= 0 )
                
                fImage(1,1)= colCorners(j,i) ;  %inverted
                fImage(1,2)= rowCorners(j,i) ;  %inverted
                fun = vgg_F_from_P(pr2, pr1) ;
                ppp2 = fImage*fun ;
                
                projMat = reshape(projMatrix(k,:,:),4,4)  ;
                
                point1 = projMat * pTemp ;
                point1 = point1 ./ point1(3,1) ;
                point1 = reshape(point1(1:2,1), 2,1) ;
                
                point2 = projMat * pTemp2 ;
                point2= point2 ./ point2(3,1) ;
                point2 = reshape(point2(1:2,1), 2,1) ;
                x1= point1(1,1) ;
                x2= point2(1,1) ;
                y1= point1(2,1) ;
                y2= point2(2,1) ;
                
                % y=ax+b
                a= (y1-y2)/(x1-x2) ;
                b= y1 - a*x1 ;
                
                color1= 0 ;
                color2=0 ;
                
                
                for x=1:1:1024
                    y = int32(((-1)*(ppp2(1,1)*x +ppp2(1,3)))/(ppp2(1,2))) ;
                    yys2(i,x) = int32(((-1)*(ppp2(1,1)*x +ppp2(1,3)))/(ppp2(1,2))) ;
                    if ( y <= 768 && y>0)
                        if ( Ib(k,y,x)  == 1 )
                            color2 = 1 ;
                        else
                            color2 = 0 ;
                        end
                        
                        if( color1 ~= color2)
                            if (color2 == 1)
                                intStarts(:,end+1) = [x ; y] ;
                                trX1 = [colCorners(j,i);rowCorners(j,i)];
                                trX2 = [x;y] ;
                                
                                XYZ1 = Triangulation(double(trX1),double(pr1), double(trX2), double(pr2)) ;
                                lastIntervals(:,end+1) = [0;coef;k;XYZ1(1,1);XYZ1(2,1);XYZ1(3,1)] ;
                                
                                
                            else
                                intEnds(:,end+1) = [x ; y] ;
                                trX1 = [colCorners(j,i);rowCorners(j,i)];
                                trX2 = [x;y] ;
                                XYZ1 = Triangulation(double(trX1),double(pr1), double(trX2), double(pr2)) ;
                                lastIntervals(:,end+1) = [0;(-1)*coef;k;XYZ1(1,1);XYZ1(2,1);XYZ1(3,1)] ;
                                
                            end
                            
                        end
                        color1 = color2 ;
                    end
                end
                
            end
            
        end
        
        lastIntervals(1,:) = sum((lastIntervals(4:6,:)- cameraCenters(1:3,lastIntervals(3,:))).^2) ;
        lastIntervals2 = intervalIntersect(lastIntervals,4) ;
        LI = reshape( lastIntervals2(4:6,:) , 3 , size(lastIntervals2,2) ) ;
        lastXw1 = lastXw ;
        lastXw(:,(end+1):(end + size(LI,2))) = LI ;
    end
    
end

figure ;
plot3((-1)*lastXw(1,:),(-1)*lastXw(2,:),(-1)*lastXw(3,:),'.') ; hold on;
plot3((-1)*cameraCenters(1,:),(-1)*cameraCenters(2,:),(-1)*cameraCenters(3,:),'+r') ;

ALL_3D = lastXw ;
ALL_3D = ALL_3D - repmat(ALL_3D(:,1),1,size(ALL_3D,2) ) ;

x = transpose([ALL_3D(1,:),0]);
y = transpose( [ALL_3D(2,:),0]);
z = transpose( [ALL_3D(3,:),0]);

DD = Delaunay3(x,y,z);
Tes = DD(:,:);
X = [x(:) y(:) z(:)];
figure ;
tetramesh(Tes,X);


