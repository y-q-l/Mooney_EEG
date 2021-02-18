
% Manipulating images
clear;

p=dir('*.jpg'); % all the file name of the images
amount=length(p); % the amount of images
picnamelist=cell(amount,2); % define new order of images
for k=1:amount
    file=fullfile(p(k).name);  % define the value of every image
    picnamelist{k,1}=p(k).name;
    picnamelist{k,2}=k;
    Img=imread(file);
    if isrgb(Img)
        Igray=rgb2gray(Img);
    else
        Igray=Img;
    end
    IgrayHis=histeq(Igray);  % histogram equalization
    [M,N]=size(IgrayHis);
    
    if (M<1000) && (N<1000) && (M>100) && (N>100) % set image size (adjustable)
        
        f=double(IgrayHis);
        g=fftshift(fft2(f));
        
        m0=fix(M/2);
        n0=fix(N/2);
        result=zeros(M,N);
        
        for Fc=10:10:80 % Fc: frequency-domain Gaussian filter (adjustable)
            for i=1:M
                for j=1:N
                    d2=(i-m0)^2+(j-n0)^2;
                    % low-pass Gaussian filter:
                    h=exp(-0.5*d2/(Fc^2));
                    result(i,j)=h*g(i,j);
                end
            end
            
            result=ifftshift(result);
            J=ifft2(result);
            X=uint8(real(J));   % X: grayscale image (low-pass)
            
            Xline=reshape(X,1,M*N);
            Xmedian=median(Xline);
            Xmediandouble=double(Xmedian);
            level= Xmediandouble/255; % median grey level
            Xbw = im2bw(X,level); % grayscale to black-white
            
            newnamelength=length(p(k).name)-4;  % to name
            name=p(k).name(1:newnamelength);
            FcStr=num2str(Fc);
            imwrite(Xbw,strcat(FcStr,'Hz','_',name, '.jpg'),'JPEG');
        end
        imwrite(IgrayHis,strcat('gray','_',name,'.jpg'),'JPEG');
    end
    clearvars -except p amount k picnamelist;
end