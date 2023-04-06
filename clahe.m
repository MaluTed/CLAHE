% function clahe() 
    filename = '/Users/malu/Desktop/红外靶标图/环境温度26度/环境温度26度/负100mK.png';
    image16bit = imread(filename);
% 
%     num_tiles_x = 4;
%     num_tiles_y = 4;
%     clip_limit = 0.1; % 尝试调整这个值
% 
%     image_double = double(image16bit) / 65535;
% 
%     adapt_histeq_img = clhe(image_double, 'NumTiles', [num_tiles_y num_tiles_x], 'ClipLimit', clip_limit, 'NBins', 65536);
%     adapt_histeq_img_16bit = uint16(adapt_histeq_img * 65535);
%     adapt_histeq_img_high8bit = uint8(bitshift(adapt_histeq_img_16bit, -8));
% 
%     figure;
%     imshow(adapt_histeq_img_high8bit, 'DisplayRange', []);
%     title(sprintf('Clip Limit: %0.2f', clip_limit));
% end

% 创建一个 256x256 的图像
% img = zeros(255, 255, 'uint16');


% % 将图像分为 3x3 网格，每块图像大小为 85x85
% tile_size = 85;
% for i = 1:3
%     for j = 1:3
%         % 为每个区域设置不同的均值，如 10000, 20000, 30000 等
%         mean_value = 10000 * (i + j - 1);
%         % 设置每个区域的正态分布
%         img((i-1)*tile_size+1:i*tile_size, (j-1)*tile_size+1:j*tile_size) = mean_value + uint16(1000 * randn(tile_size, tile_size));
%     end
% end

img = image16bit;

% ----------------------
% 显示原始图像
figure;
subplot(4,3,1);
imshow(img, 'DisplayRange', [min(img(:)) max(img(:))], 'InitialMagnification', 'fit');
title('16bit原始图像，7x7,bins=1024');

img_normalized = double(img) / double(max(img(:)));
% 显示原始图像的3D图
subplot(4,3,2);
mesh(img_normalized);
title('原始图像的3D图');

num_tiles_x = 4;
num_tiles_y = 4;
clip_limit = 0.5;




subplot(4,3,3);
imhist(img_normalized);
title('原始图像的直方图');





% ----------------------
% 使用 histeq 处理图像
img_histeq = histeq(img);

% 显示处理后的图像
subplot(4,3,4);
imshow(img_histeq, 'DisplayRange', [min(img_histeq(:)) max(img_histeq(:))], 'InitialMagnification', 'fit');
title('经过 histeq 处理的图像');


img_histeq_normalized = double(img_histeq) / double(max(img_histeq(:)));
subplot(4,3,5);
mesh(img_histeq)
title('经过 histeq 处理的3D图')


subplot(4,3,6);
imhist(img_histeq_normalized);
title('经过 histeq 处理的直方图')

% ----------------------
% 显示处理后的图像
% 使用 adapthisteq 处理图像
img_adapthisteq = adapthisteq(img, 'NumTiles', [num_tiles_y num_tiles_x], 'ClipLimit', clip_limit, 'NBins', 65536);
subplot(4,3,7);
imshow(img_adapthisteq, 'DisplayRange', [min(img_adapthisteq(:)) max(img_adapthisteq(:))], 'InitialMagnification', 'fit');
title('经过 adapthisteq 处理的图像');

% 为了显示处理后的3D图和直方图，我们需要将图像转换为 double 类型并归一化到 0-1 范围内
img_adapthisteq_normalized = double(img_adapthisteq) / double(max(img_adapthisteq(:)));

% 显示处理后的3D图
subplot(4,3,8);
mesh(img_adapthisteq_normalized);
title('经过 adapthisteq 处理的3D图');

% 显示处理后的直方图
subplot(4,3,9);
imhist(img_adapthisteq_normalized);
title('经过 adapthisteq 处理的的直方图');

% ----------------------

img_adapthisteq = clhe(img, 'NumTiles', [num_tiles_y num_tiles_x], 'ClipLimit', clip_limit, 'NBins', 65536);
subplot(4,3,10);
imshow(img_adapthisteq, 'DisplayRange', [min(img_adapthisteq(:)) max(img_adapthisteq(:))], 'InitialMagnification', 'fit');
title('经过 clhe 处理的图像');

img_adapthisteq_normalized = double(img_adapthisteq) / double(max(img_adapthisteq(:)));

subplot(4,3,11);
mesh(img_adapthisteq_normalized);
title('经过 clhe 处理的3D图');

subplot(4,3,12);
imhist(img_adapthisteq_normalized);
title('经过 clhe 处理的的直方图');

