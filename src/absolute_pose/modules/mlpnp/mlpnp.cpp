/******************************************************************************
* Author:   Steffen Urban                                              *
* Contact:  urbste@gmail.com                                          *
* License:  Copyright (c) 2016 Steffen Urban, ANU. All rights reserved.      *
*                                                                            *
* Redistribution and use in source and binary forms, with or without         *
* modification, are permitted provided that the following conditions         *
* are met:                                                                   *
* * Redistributions of source code must retain the above copyright           *
*   notice, this list of conditions and the following disclaimer.            *
* * Redistributions in binary form must reproduce the above copyright        *
*   notice, this list of conditions and the following disclaimer in the      *
*   documentation and/or other materials provided with the distribution.     *
* * Neither the name of ANU nor the names of its contributors may be         *
*   used to endorse or promote products derived from this software without   *
*   specific prior written permission.                                       *
*                                                                            *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
* ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
* SUCH DAMAGE.                                                               *
******************************************************************************/

#include <opengv/absolute_pose/modules/mlpnp.hpp>

void
opengv::absolute_pose::modules::mlpnp::mlpnpJacs(
	const point_t& pt,
	const Eigen::Vector3d& nullspace_r,
	const Eigen::Vector3d& nullspace_s,
	const cayley_t& c,
	const translation_t& t,
Eigen::MatrixXd& jacs)
{
	double r1 = nullspace_r[0];
	double r2 = nullspace_r[1];
	double r3 = nullspace_r[2];

	double s1 = nullspace_s[0];
	double s2 = nullspace_s[1];
	double s3 = nullspace_s[2];

	double X1 = pt[0];
	double X2 = pt[1];
	double X3 = pt[2];

	double c1 = c[0];
	double c2 = c[1];
	double c3 = c[2];

	double T1 = t[0];
	double T2 = t[1];
	double T3 = t[2];

	double	t2 = c1*c1;
	double	t3 = c2*c2;
	double	t4 = c3*c3;
	double	t10 = T2*c3*2.0;
	double	t11 = T3*c2*2.0;
	double	t12 = X2*c3*2.0;
	double	t13 = X3*c2*2.0;
	double	t14 = T1*t2;
	double	t15 = T1*t3;
	double	t16 = T1*t4;
	double	t17 = X1*t2;
	double	t18 = X1*t3;
	double	t19 = X1*t4;
	double	t20 = T2*c1*c2*2.0;
	double	t21 = T3*c1*c3*2.0;
	double	t22 = X2*c1*c2*2.0;
	double	t23 = X3*c1*c3*2.0;
	double	t5 = T1 - X1 + t10 - t11 - t12 + t13 + t14 - t15 - t16 - t17 + t18 + t19 + t20 + t21 - t22 - t23;
	double	t6 = t2 + t3 + t4 + 1.0;
	double	t7 = 1.0 / (t6*t6);
	double	t26 = T1*c3*2.0;
	double	t27 = T3*c1*2.0;
	double	t28 = X1*c3*2.0;
	double	t29 = X3*c1*2.0;
	double	t30 = T2*t2;
	double	t31 = T2*t3;
	double	t32 = T2*t4;
	double	t33 = X2*t2;
	double	t34 = X2*t3;
	double	t35 = X2*t4;
	double	t36 = T1*c1*c2*2.0;
	double	t37 = T3*c2*c3*2.0;
	double	t38 = X1*c1*c2*2.0;
	double	t39 = X3*c2*c3*2.0;
	double	t8 = T2 - X2 - t26 + t27 + t28 - t29 - t30 + t31 - t32 + t33 - t34 + t35 + t36 + t37 - t38 - t39;
	double	t42 = T1*c2*2.0;
	double	t43 = T2*c1*2.0;
	double	t44 = X1*c2*2.0;
	double	t45 = X2*c1*2.0;
	double	t46 = T3*t2;
	double	t47 = T3*t3;
	double	t48 = T3*t4;
	double	t49 = X3*t2;
	double	t50 = X3*t3;
	double	t51 = X3*t4;
	double	t52 = T1*c1*c3*2.0;
	double	t53 = T2*c2*c3*2.0;
	double	t54 = X1*c1*c3*2.0;
	double	t55 = X2*c2*c3*2.0;
	double	t9 = T3 - X3 + t42 - t43 - t44 + t45 - t46 - t47 + t48 + t49 + t50 - t51 + t52 + t53 - t54 - t55;
	double	t24 = t5*t5;
	double	t25 = t7*t24;
	double	t40 = t8*t8;
	double	t41 = t7*t40;
	double	t56 = t9*t9;
	double	t57 = t7*t56;
	double	t58 = t25 + t41 + t57;
	double	t59 = 1.0 / sqrt(t58);
	double	t60 = 1.0 / t6;
	double	t61 = T3*2.0;
	double	t62 = X3*2.0;
	double	t63 = t42 - t43 - t44 + t45 + t61 - t62;
	double	t64 = T2*2.0;
	double	t65 = X2*2.0;
	double	t66 = t26 - t27 - t28 + t29 - t64 + t65;
	double	t67 = 1.0 / (t6*t6*t6);
	double	t68 = T1*c1*2.0;
	double	t69 = T2*c2*2.0;
	double	t70 = T3*c3*2.0;
	double	t76 = X1*c1*2.0;
	double	t77 = X2*c2*2.0;
	double	t78 = X3*c3*2.0;
	double	t71 = t68 + t69 + t70 - t76 - t77 - t78;
	double	t72 = 1.0 / pow(t58, 3.0 / 2.0);
	double	t73 = c1*t24*t67*4.0;
	double	t74 = c1*t40*t67*4.0;
	double	t75 = c1*t56*t67*4.0;
	double	t80 = t7*t8*t63*2.0;
	double	t81 = t7*t9*t66*2.0;
	double	t82 = t5*t7*t71*2.0;
	double	t79 = t73 + t74 + t75 - t80 - t81 - t82;
	double	t83 = T1*2.0;
	double	t84 = X1*2.0;
	double	t85 = t10 - t11 - t12 + t13 + t83 - t84;
	double	t86 = t5*t7*t63*2.0;
	double	t87 = c2*t24*t67*4.0;
	double	t88 = c2*t40*t67*4.0;
	double	t89 = c2*t56*t67*4.0;
	double	t91 = t7*t9*t85*2.0;
	double	t92 = t7*t8*t71*2.0;
	double	t90 = t86 + t87 + t88 + t89 - t91 - t92;
	double	t93 = t7*t8*t85*2.0;
	double	t94 = c3*t24*t67*4.0;
	double	t95 = c3*t40*t67*4.0;
	double	t96 = c3*t56*t67*4.0;
	double	t97 = t5*t7*t66*2.0;
	double	t99 = t7*t9*t71*2.0;
	double	t98 = t93 + t94 + t95 + t96 + t97 - t99;
	double	t100 = c3*2.0;
	double	t106 = c1*c2*2.0;
	double	t101 = t100 - t106;
	double	t102 = c2*2.0;
	double	t103 = c1*c3*2.0;
	double	t104 = t102 + t103;
	double	t105 = t2 - t3 - t4 + 1.0;
	double	t107 = t7*t9*t104*2.0;
	double	t108 = t5*t7*t105*2.0;
	double	t110 = t7*t8*t101*2.0;
	double	t109 = t107 + t108 - t110;
	double	t111 = t100 + t106;
	double	t112 = c1*2.0;
	double	t115 = c2*c3*2.0;
	double	t113 = t112 - t115;
	double	t114 = t2 - t3 + t4 - 1.0;
	double	t116 = t7*t9*t113*2.0;
	double	t117 = t7*t8*t114*2.0;
	double	t119 = t5*t7*t111*2.0;
	double	t118 = t116 + t117 - t119;
	double	t120 = t102 - t103;
	double	t121 = t112 + t115;
	double	t122 = t2 + t3 - t4 - 1.0;
	double	t123 = t5*t7*t120*2.0;
	double	t124 = t7*t9*t122*2.0;
	double	t126 = t7*t8*t121*2.0;
	double	t125 = t123 + t124 - t126;

	jacs(0, 0) = -r2*t59*t60*t63 - r3*t59*t60*t66 - r1*t59*t60*t71 + c1*r1*t5*t7*t59*2.0 + c1*r2*t7*t8*t59*2.0 + c1*r3*t7*t9*t59*2.0 - r1*t5*t60*t72*t79*(1.0 / 2.0) - r2*t8*t60*t72*t79*(1.0 / 2.0) - r3*t9*t60*t72*t79*(1.0 / 2.0);
	jacs(0, 1) = r1*t59*t60*t63 - r2*t59*t60*t71 - r3*t59*t60*t85 + c2*r1*t5*t7*t59*2.0 + c2*r2*t7*t8*t59*2.0 + c2*r3*t7*t9*t59*2.0 - r1*t5*t60*t72*t90*(1.0 / 2.0) - r2*t8*t60*t72*t90*(1.0 / 2.0) - r3*t9*t60*t72*t90*(1.0 / 2.0);
	jacs(0, 2) = -r3*t59*t60*t71 + r2*t59*t60*t85 + r1*t59*t60*(t26 - t27 - t28 + t29 - t64 + t65) + c3*r1*t5*t7*t59*2.0 + c3*r2*t7*t8*t59*2.0 + c3*r3*t7*t9*t59*2.0 - r1*t5*t60*t72*t98*(1.0 / 2.0) - r2*t8*t60*t72*t98*(1.0 / 2.0) - r3*t9*t60*t72*t98*(1.0 / 2.0);
	jacs(0, 3) = r2*t59*t60*t101 - r1*t59*t60*t105 - r3*t59*t60*t104 + r1*t5*t60*t72*t109*(1.0 / 2.0) + r2*t8*t60*t72*t109*(1.0 / 2.0) + r3*t9*t60*t72*t109*(1.0 / 2.0);
	jacs(0, 4) = -r1*t59*t60*t111 + r2*t59*t60*t114 + r3*t59*t60*t113 - r1*t5*t60*t72*t118*(1.0 / 2.0) - r2*t8*t60*t72*t118*(1.0 / 2.0) - r3*t9*t60*t72*t118*(1.0 / 2.0);
	jacs(0, 5) = r1*t59*t60*t120 - r2*t59*t60*t121 + r3*t59*t60*t122 - r1*t5*t60*t72*t125*(1.0 / 2.0) - r2*t8*t60*t72*t125*(1.0 / 2.0) - r3*t9*t60*t72*t125*(1.0 / 2.0);
	jacs(1, 0) = -s2*t59*t60*t63 - s3*t59*t60*t66 - s1*t59*t60*t71 + c1*s1*t5*t7*t59*2.0 + c1*s2*t7*t8*t59*2.0 + c1*s3*t7*t9*t59*2.0 - s1*t5*t60*t72*t79*(1.0 / 2.0) - s2*t8*t60*t72*t79*(1.0 / 2.0) - s3*t9*t60*t72*t79*(1.0 / 2.0);
	jacs(1, 1) = s1*t59*t60*t63 - s2*t59*t60*t71 - s3*t59*t60*t85 + c2*s1*t5*t7*t59*2.0 + c2*s2*t7*t8*t59*2.0 + c2*s3*t7*t9*t59*2.0 - s1*t5*t60*t72*t90*(1.0 / 2.0) - s2*t8*t60*t72*t90*(1.0 / 2.0) - s3*t9*t60*t72*t90*(1.0 / 2.0);
	jacs(1, 2) = -s3*t59*t60*t71 + s2*t59*t60*t85 + s1*t59*t60*(t26 - t27 - t28 + t29 - t64 + t65) + c3*s1*t5*t7*t59*2.0 + c3*s2*t7*t8*t59*2.0 + c3*s3*t7*t9*t59*2.0 - s1*t5*t60*t72*t98*(1.0 / 2.0) - s2*t8*t60*t72*t98*(1.0 / 2.0) - s3*t9*t60*t72*t98*(1.0 / 2.0);
	jacs(1, 3) = s2*t59*t60*t101 - s1*t59*t60*t105 - s3*t59*t60*t104 + s1*t5*t60*t72*t109*(1.0 / 2.0) + s2*t8*t60*t72*t109*(1.0 / 2.0) + s3*t9*t60*t72*t109*(1.0 / 2.0);
	jacs(1, 4) = -s1*t59*t60*t111 + s2*t59*t60*t114 + s3*t59*t60*t113 - s1*t5*t60*t72*t118*(1.0 / 2.0) - s2*t8*t60*t72*t118*(1.0 / 2.0) - s3*t9*t60*t72*t118*(1.0 / 2.0);
	jacs(1, 5) = s1*t59*t60*t120 - s2*t59*t60*t121 + s3*t59*t60*t122 - s1*t5*t60*t72*t125*(1.0 / 2.0) - s2*t8*t60*t72*t125*(1.0 / 2.0) - s3*t9*t60*t72*t125*(1.0 / 2.0);
}