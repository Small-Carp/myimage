from PIL import Image
import numpy as np
# import scipy
import matplotlib.pyplot as plt
import math

target_size = (400, 400)
flag = np.zeros((200, 200))
# 边界距离矩阵
flag_border = np.zeros((200, 200))
flag_next1 = np.zeros((200, 200))
Threshold = 30 # 阈值
emptyImage = np.zeros(target_size)  # 新的图
th, tw = target_size[0], target_size[1]

def search(data,x1,y1,x2,y2):
    for r in range(x1-1,x1+2):
        for w in range(y1 - 1, y1 + 2):
            if(r<200 and  w<200):
                if((r!=x1 and w!=y1)and(r!=x2 and w!=y2) and flag_border[r,w]==1):
                    return data[r,w]

    return (data[x1,y1]+data[x2,y2])*1.0/2


def Linear_compensation(data,i, j, real_x,real_y,x1,y1,x2,y2):#边界补偿
    k = math.hypot(abs(x1 - x2), abs(y1 - y2))
    v = math.hypot(abs(real_x - x1), abs(real_y - y1))
    #寻找x1,y1,x2,y2的邻近点
    datx1=search(data,x1,y1,x2,y2)
    datx2 = search(data, x2, y2, x1, y1)
    RL=data[x1,y1]-(datx1+data[x2,y2])/2
    RR=data[x2,y2]-(data[x1,y1]+datx2)/2
    R1VK=4*(data[x1,y1]-(datx1+data[x2,y2])/2)*v/k*(1-v/k)
    R2VK=4*(data[x2,y2]-(data[x1,y1]+datx2)/2)*v/k*(1-v/k)
    if(RL*RR>0 and abs(RL)<=abs(RR)):
        modR=RL
    elif(RL*RR>0 and abs(RL)>=abs(RR)):
        modR=RR
    else:
        modR=0
    if(modR==RL):
        MIJ=R1VK
    elif(modR==RR):
        MIJ=R2VK
    else:
        MIJ=0

    #t1 = math.hypot(abs(i - p1x), abs(j - p1y))

    partvalue=(1.0-v*1.0/k)*data[x1,y1]+(v*1.0/k)*data[x2,y2]

    compensatevalue=0.25*MIJ

    emptyImage[i, j]=partvalue+compensatevalue



# 一点放大

def one_point(data,i, j, real_x,real_y,point1,point2,point3,point4):
    if (flag_border[point1[0], point1[1]]==0):
        emptyImage[i,j]=data[point1[0], point1[1]]
    if (flag_border[point2[0], point2[1]]==0):
        emptyImage[i, j] = data[point2[0], point2[1]]
    if (flag_border[point3[0], point3[1]]==0):
        emptyImage[i, j] = data[point3[0], point3[1]]
    if (flag_border[point4[0], point4[1]]==0):
        emptyImage[i, j] = data[point4[0], point4[1]]





def two_point(data,i, j, real_x,real_y,point1,point2,point3,point4):
    zs,ys,zx,yx=0,0,0,0
    tem = 0
    if (flag_border[point1[0], point1[1]] == 1):
        zs = 1
        tem = tem + 1
    if (flag_border[point2[0], point2[1]] == 1):
        ys = 1
        tem = tem + 1
    if (flag_border[point3[0], point3[1]] == 1):
        zx = 1
        tem = tem + 1
    if (flag_border[point4[0], point4[1]] == 1):
        yx = 1
        tem = tem + 1
    ans1 = 0
    if (zs == 1 and ys == 1):
        v = math.hypot(abs(real_x - point3[0]), abs(real_y - point3[1]))
        u = math.hypot(abs(real_x - point4[0]), abs(real_y - point4[1]))
        tem1=v/(v+u)*data[point4[0],point4[1]]
        tem2 = u / (v + u) * data[point3[0], point3[1]]
        emptyImage[i, j] = tem1+tem2
    if (zs == 1 and zx == 1):
        v = math.hypot(abs(real_x - point2[0]), abs(real_y - point2[1]))
        u = math.hypot(abs(real_x - point4[0]), abs(real_y - point4[1]))
        tem1 = v / (v + u) * data[point4[0], point4[1]]
        tem2 = u / (v + u) * data[point2[0], point2[1]]
        emptyImage[i, j] = tem1 + tem2
    if (zx == 1 and yx == 1):
        v = math.hypot(abs(real_x - point1[0]), abs(real_y - point1[1]))
        u = math.hypot(abs(real_x - point2[0]), abs(real_y - point2[1]))
        tem1 = v / (v + u) * data[point2[0], point2[1]]
        tem2 = u / (v + u) * data[point1[0], point1[1]]
        emptyImage[i, j] = tem1 + tem2
    if (ys == 1 and yx== 1):
        v = math.hypot(abs(real_x - point1[0]), abs(real_y - point1[1]))
        u = math.hypot(abs(real_x - point3[0]), abs(real_y - point3[1]))
        tem1 = v / (v + u) * data[point3[0], point3[1]]
        tem2 = u / (v + u) * data[point1[0], point1[1]]
        emptyImage[i, j] = tem1 + tem2



# 三点放大
def three_point(data,i, j, real_x,real_y,point1,point2,point3,point4):
    if(flag_border[point1[0],point1[1]]):
        p=(point4[0]-real_x)*1.0
        q=(point4[1]-real_y)*1.0
        temp1=p*data[point3[0],point3[1]]
        temp2 = q* data[point2[0], point2[1]]
        temp3 = (1-p-q) * data[point4[0], point4[1]]
        emptyImage[i,j]=temp1+temp2+temp3
    if (flag_border[point2[0], point2[1]]):
        p = (real_x-point3[0]) * 1.0
        q = (point3[1] - real_y) * 1.0
        temp1 = p * data[point4[0], point4[1]]
        temp2 = q * data[point1[0], point1[1]]
        temp3 = (1 - p - q) * data[point3[0], point3[1]]
        emptyImage[i, j] = temp1 + temp2 + temp3
    if (flag_border[point3[0], point3[1]]):
        p = (point2[0] - real_x) * 1.0
        q = (real_y-point2[1]  ) * 1.0
        temp1 = p * data[point1[0], point1[1]]
        temp2 = q * data[point4[0], point4[1]]
        temp3 = (1 - p - q) * data[point2[0], point2[1]]
        emptyImage[i, j] = temp1 + temp2 + temp3

    if (flag_border[point4[0], point4[1]]):
        p = (real_x-point1[0] ) * 1.0
        q = (real_y-point1[1]) * 1.0
        temp1 = p * data[point2[0], point2[1]]
        temp2 = q * data[point3[0], point3[1]]
        temp3 = (1 - p - q) * data[point1[0], point1[1]]
        emptyImage[i, j] = temp1 + temp2 + temp3


# 四点放大
def bi_linear(data,i, j, real_x,real_y,point1,point2,point3,point4):
    # 读取输入图像
    x=float(real_x-point1[0])
    y=float(real_y-point1[1])
    t1=(1-x)*(1-y)*data[point1[0],point1[1]]
    t2= (1 - x)*y * data[point2[0], point2[1]]
    t3 = x*(1-y)* data[point3[0], point3[1]]
    t4 = x*y * data[point4[0], point4[1]]
    emptyImage[i,j]=t1+t2+t3+t4




#     emptyImage.show()
#     emptyImage.save( 'E:/cat.jpg')
# 边界
# 边界程度
# 边界距离
# 补偿项

# def deal():
def Boundary_extraction(data):  # 计算边界距离
    # 获得矩阵行列信息
    row_len = data.shape[0]#行
    Column_len = data.shape[1]#列
    # print( row_len,Column_len,data[1][1],data[2][1])
    for r in range(row_len):
        for c in range(Column_len):
            if r == 0 or c == 0 or r == row_len - 1 or c == Column_len - 1:#排除图像边界影响
                flag[r, c] = 0

            else:
                a = abs(data[r, c] - data[r - 1, c])  # 上
                b = abs(data[r, c] - data[r + 1, c])  # 下
                c = abs(data[r, c] - data[r, c - 1])  # 左

                d = abs(data[r, c] - data[r, c + 1])  # 右
                flag[r, c] = max(a, b, c, d)#得到边界的最大绝对值
                if flag[r, c]*0.2 > Threshold:
                    flag_border[r, c] = 1  # 此点为边界点

    print(flag)
    print(flag_border)

    # print(row_len,Column_len)


# mian()
def ImageToMatrix(src):#得到灰度矩阵
    # 读取图片
    im = Image.open(src)
    # 显示图片
    im.show()
    width, height = im.size
    im = im.convert("L")#使用L方式打开为灰度图像
    data = im.getdata()  # 灰度图像
    # data = np.matrix(data,dtype='float')/255.0
    #data = np.matrix(data, dtype='int')#转化为整型灰度矩阵
    data = np.matrix(data, dtype='int')
    # new_data = np.reshape(data,(width,height))
    new_data = np.reshape(data, (height, width))
    return new_data


#     new_im = Image.fromarray(new_data)
#     # 显示图片
#     new_im.show()
def MatrixToImage(data):
    # data = data*255
    new_im = Image.fromarray(data.astype(np.uint8))
    return new_im

def analysis(i,j,p1x,p1y,p2x,p2y):
    t1=math.hypot(abs(i-p1x),abs(j-p1y))
    t2 = math.hypot(abs(i - p2x), abs(j - p2y))
    t3=t1+t2
    t4=math.hypot(abs(p1x - p2x), abs(p1y- p2y))
    if(t3==t4):
        return 1
    else:
        return 0


# 综合处理
def amplification(data):
    row_len = data.shape[0]
    Column_len = data.shape[1]
    for i in range(th):#th,tw预先定义了
        for j in range(tw):

            #             corr_x = (i+0.5)/th*pic.shape[0]-0.5
            #             corr_y = (j+0.5)/tw*pic.shape[1]-0.5
            real_x = float((i+0.0) / th * row_len+0.0)
            real_y = float((j+0.0) / tw * Column_len+0.0)
            # if i*pic.shape[0]%th==0 and j*pic.shape[1]%tw==0:     # 对应的点正好是一个像素点，直接拷贝
            #   emptyImage[i, j, k] = pic[int(corr_x), int(corr_y), k]
            point1 = (math.floor(real_x), math.floor(real_y))  # 左上角的点
            point2 = (point1[0], point1[1] + 1)  # 右上角点
            point3 = (point1[0] + 1, point1[1])  # 左下角
            point4 = (point1[0] + 1, point1[1] + 1)  # 右下角
            if(point1[0]<200 and point1[1]<200 and point2[0]<200 and point2[1]<200 and point3[0]<200 and point3[1]<200 and point4[0]<200 and point4[1]<200  ):
                if(real_x==point1[0] and real_y==point1[1]):
                    emptyImage[i, j] = data[point1[0], point1[1]]
                else:
                    tem=0
                    zs,ys,zx,yx=0,0,0,0
                    if(flag_border[point1[0],point1[1]] == 1):
                        zs=1
                        tem=tem+1
                    if (flag_border[point2[0], point2[1]] == 1):
                        ys = 1
                        tem = tem + 1
                    if (flag_border[point3[0], point3[1]] == 1):
                        zx = 1
                        tem = tem + 1
                    if (flag_border[point4[0], point4[1]] == 1):
                        yx = 1
                        tem = tem + 1
                    if(tem==0):#远离边界的平坦区
                        bi_linear(data,i, j, real_x,real_y,point1,point2,point3,point4)
                    if(tem==1):#衔接
                        #三角形
                        three_point(data,i, j, real_x,real_y,point1,point2,point3,point4)
                    if (tem == 2):
                        #判断是否在边缘线上
                        ans1 = 0
                        if (zs == 1 and ys == 1):#1,2点
                            ans1 = analysis(real_x,real_y, point1[0], point1[1], point2[0], point2[1]) + ans1
                            if (ans1):
                                Linear_compensation(data, i, j, real_x,real_y,point1[0], point1[1], point2[0], point2[1])
                        if (zs == 1 and zx == 1):#1,3点
                            ans1 = analysis(real_x, real_y, point1[0], point1[1], point3[0], point3[1]) + ans1
                            if (ans1):
                                Linear_compensation(data, i, j, real_x,real_y,point1[0], point1[1], point3[0], point3[1])
                        if (zs == 1 and yx == 1):#1,4点
                            ans1 = analysis(real_x, real_y, point1[0], point1[1], point4[0], point4[1]) + ans1
                            if (ans1):
                                Linear_compensation(data, i, j, real_x,real_y,point1[0], point1[1], point4[0], point4[1])
                        if (zx == 1 and ys == 1):#2，3点
                            ans1 = analysis(real_x, real_y, point2[0], point2[1], point3[0], point3[1]) + ans1
                            if (ans1):
                                Linear_compensation(data, i, j,real_x,real_y, point2[0], point2[1], point3[0], point3[1])
                        if (yx == 1 and ys == 1):#2,4点
                            ans1 = analysis(real_x, real_y, point2[0], point2[1], point4[0], point4[1]) + ans1
                            if (ans1):
                                Linear_compensation(data, i, j,real_x,real_y, point2[0], point2[1], point4[0], point4[1])
                        if (zx == 1 and yx == 1):#3,4点
                            ans1 = analysis(real_x, real_y, point3[0], point3[1], point4[0], point4[1]) + ans1
                            if (ans1):
                                Linear_compensation(data, i, j,real_x,real_y, point3[0], point3[1], point4[0], point4[1])

                        if(ans1==0):
                        # 三角形
                            two_point(data,i, j, real_x,real_y,point1,point2,point3,point4)

                    if (tem == 3):
                            # 最近邻
                        ans1=0
                        if (zs == 1 and ys == 1):#1,2
                            ans1 = analysis(real_x, real_y, point1[0], point1[1], point2[0], point2[1]) + ans1
                            if(ans1):
                                Linear_compensation(data, i, j,real_x,real_y, point1[0], point1[1], point2[0], point2[1])
                        if (zs == 1 and zx == 1):#1,3
                            ans1 = analysis(real_x, real_y, point1[0], point1[1], point3[0], point3[1]) + ans1
                            if (ans1):
                                Linear_compensation(data,i, j, real_x,real_y,point1[0], point1[1], point3[0], point3[1])
                        if (zs == 1 and yx == 1):#1,4
                            ans1 = analysis(real_x, real_y, point1[0], point1[1], point4[0], point4[1]) + ans1
                            if (ans1):
                                Linear_compensation(data,i, j,real_x,real_y, point1[0], point1[1], point4[0], point4[1])
                        if (ys == 1 and zx == 1):#2,3
                            ans1 = analysis(real_x, real_y, point2[0], point3[1], point3[0], point3[1]) + ans1
                            if (ans1):
                                Linear_compensation(data,i, j,real_x,real_y, point2[0], point2[1], point3[0], point3[1])
                        if (ys== 1 and yx== 1):#2,4
                            ans1 = analysis(real_x, real_y, point2[0], point2[1], point4[0], point4[1]) + ans1
                            if (ans1):
                                Linear_compensation(data,i, j,real_x,real_y, point2[0], point2[1], point4[0], point4[1])
                        if (zx == 1 and yx == 1):#3,4
                            ans1 = analysis(real_x, real_y, point3[0], point3[1], point4[0], point4[1]) + ans1
                            if (ans1):
                                Linear_compensation(data,i, j,real_x,real_y, point3[0], point3[1], point4[0], point4[1])
                        if ans1==0:
                            one_point(data,i, j, real_x,real_y,point1,point2,point3,point4)
    return data





def main():
    src = 'C:/Users/chang/Desktop/cat.jpg'
    dst = 'C:/Users/chang/Desktop/cutecat.jpg'
    data = ImageToMatrix(src)  # 灰度矩阵第一次显示原图

    old_im=MatrixToImage(data)
    #plt.imshow(data, cmap=plt.cm.gray, interpolation='nearest')
    old_im.show()#旧灰度图
    print(data)

    Boundary_extraction(data)
    deal_data=amplification(data)#处理完data矩阵
    print(deal_data.shape[0])
    target_size = (400, 400)
    # target_size = (200, 200, 3)
    # def deal(src,dst,target_size)
    #new_im = MatrixToImage(data)
    new_im = MatrixToImage(emptyImage)
    #plt.imshow(data, cmap=plt.cm.gray, interpolation='nearest')
    #plt.imshow(deal_data, cmap=plt.cm.gray, interpolation='nearest')
    new_im.show()#新灰度图
    new_im.save(dst)
    # data1 = bi_linear(data, target_size)
    # new_im1 = MatrixToImage(data1)
    # # plt.imshow(data1, cmap=plt.cm.gray, interpolation='nearest')
    # new_im1.show()


if __name__ == '__main__':
    main()