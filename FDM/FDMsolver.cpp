#ifndef _FDMSOLVER_HPP_
#define _FDMSOLVER_HPP_

#include "FDMsolver.h"

void FDMsolver::Calculator(){
	// 1. 按参数nx, ny, nz划分宏模型表面，形成2(nx*ny + ny*nz + nz*nx)块边界面元
	// 2. 取定每个面元中心点，连接两两相对的面元中心点，形成立体网格
	// 3. 对内部节点根据拉普拉斯方程列式，对面元中心点根据dU/dx = -E导函数关系列式，整理得到分块矩阵
	// 4. 如果含有内部导体，则需要对矩阵进行压缩
	// 5. 通过对分块矩阵各个分块进行矩阵计算，得到消去内部节点的电容矩阵
	// 6. 打印返回 
}

#endif // _FDMSOLVER_HPP_