# 报告

用各种方法对 Runge 函数插值, 结果如下. 

![interpolations](interpolations.jpg)

可以看出: 均匀取点的牛顿法有很大的振荡偏差, 对于 Runge 函数效果很差; 线性分段插值虽然没有很大的偏差, 但是精度不足; 不均匀取点的拉格朗日法以及均匀取点的分段立方 Hermite 插值效果则很好.