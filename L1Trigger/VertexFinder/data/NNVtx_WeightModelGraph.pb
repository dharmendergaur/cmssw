
@
weightPlaceholder*
dtype0*
shape:���������
K
model_9/q_activation_8/mul_3/xConst*
dtype0*
valueB
 *  �?
K
model_9/q_activation_7/mul_3/xConst*
dtype0*
valueB
 *  �?
K
model_9/q_activation_6/mul_3/xConst*
dtype0*
valueB
 *  �?
�
(model_9/weight_1/ReadVariableOp/resourceConst*
dtype0*�
value�B�
"x!-C?   �   ��g��    ����   �   ���j����/�ɢ/?�D�       ��   �    di�>;s�>ҏ>>>��>���>    "�?           �   �
`
!model_9/weight_1/ReadVariableOp_2Identity(model_9/weight_1/ReadVariableOp/resource*
T0
E
model_9/weight_1/mul_3/xConst*
dtype0*
valueB
 *  �?
`
!model_9/weight_1/ReadVariableOp_1Identity(model_9/weight_1/ReadVariableOp/resource*
T0
I
model_9/weight_1/Neg_1Neg!model_9/weight_1/ReadVariableOp_1*
T0
E
model_9/weight_1/mul_2/xConst*
dtype0*
valueB
 *  �?
^
model_9/weight_1/ReadVariableOpIdentity(model_9/weight_1/ReadVariableOp/resource*
T0
C
model_9/weight_1/mul/yConst*
dtype0*
valueB
 *   D
]
model_9/weight_1/mulMulmodel_9/weight_1/ReadVariableOpmodel_9/weight_1/mul/y*
T0
@
model_9/weight_1/Pow/xConst*
dtype0*
value	B :
@
model_9/weight_1/Pow/yConst*
dtype0*
value	B :
T
model_9/weight_1/PowPowmodel_9/weight_1/Pow/xmodel_9/weight_1/Pow/y*
T0
[
model_9/weight_1/CastCastmodel_9/weight_1/Pow*

DstT0*

SrcT0*
Truncate( 
Y
model_9/weight_1/truedivRealDivmodel_9/weight_1/mulmodel_9/weight_1/Cast*
T0
>
model_9/weight_1/NegNegmodel_9/weight_1/truediv*
T0
B
model_9/weight_1/RoundRoundmodel_9/weight_1/truediv*
T0
T
model_9/weight_1/addAddV2model_9/weight_1/Negmodel_9/weight_1/Round*
T0
L
model_9/weight_1/StopGradientStopGradientmodel_9/weight_1/add*
T0
a
model_9/weight_1/add_1AddV2model_9/weight_1/truedivmodel_9/weight_1/StopGradient*
T0
U
(model_9/weight_1/clip_by_value/Minimum/yConst*
dtype0*
valueB
 * ��C
|
&model_9/weight_1/clip_by_value/MinimumMinimummodel_9/weight_1/add_1(model_9/weight_1/clip_by_value/Minimum/y*
T0
M
 model_9/weight_1/clip_by_value/yConst*
dtype0*
valueB
 *   �
|
model_9/weight_1/clip_by_valueMaximum&model_9/weight_1/clip_by_value/Minimum model_9/weight_1/clip_by_value/y*
T0
]
model_9/weight_1/mul_1Mulmodel_9/weight_1/Castmodel_9/weight_1/clip_by_value*
T0
I
model_9/weight_1/truediv_1/yConst*
dtype0*
valueB
 *   D
d
model_9/weight_1/truediv_1RealDivmodel_9/weight_1/mul_1model_9/weight_1/truediv_1/y*
T0
\
model_9/weight_1/mul_2Mulmodel_9/weight_1/mul_2/xmodel_9/weight_1/truediv_1*
T0
X
model_9/weight_1/add_2AddV2model_9/weight_1/Neg_1model_9/weight_1/mul_2*
T0
X
model_9/weight_1/mul_3Mulmodel_9/weight_1/mul_3/xmodel_9/weight_1/add_2*
T0
P
model_9/weight_1/StopGradient_1StopGradientmodel_9/weight_1/mul_3*
T0
l
model_9/weight_1/add_3AddV2!model_9/weight_1/ReadVariableOp_2model_9/weight_1/StopGradient_1*
T0
p
model_9/weight_1/MatMulMatMulweightmodel_9/weight_1/add_3*
T0*
transpose_a( *
transpose_b( 

*model_9/weight_1/ReadVariableOp_3/resourceConst*
dtype0*=
value4B2
"(�l��j�*Ӽ���    ���        ��<PI<
b
!model_9/weight_1/ReadVariableOp_5Identity*model_9/weight_1/ReadVariableOp_3/resource*
T0
E
model_9/weight_1/mul_7/xConst*
dtype0*
valueB
 *  �?
b
!model_9/weight_1/ReadVariableOp_4Identity*model_9/weight_1/ReadVariableOp_3/resource*
T0
I
model_9/weight_1/Neg_3Neg!model_9/weight_1/ReadVariableOp_4*
T0
E
model_9/weight_1/mul_6/xConst*
dtype0*
valueB
 *  �?
b
!model_9/weight_1/ReadVariableOp_3Identity*model_9/weight_1/ReadVariableOp_3/resource*
T0
E
model_9/weight_1/mul_4/yConst*
dtype0*
valueB
 *   D
c
model_9/weight_1/mul_4Mul!model_9/weight_1/ReadVariableOp_3model_9/weight_1/mul_4/y*
T0
B
model_9/weight_1/Pow_1/xConst*
dtype0*
value	B :
B
model_9/weight_1/Pow_1/yConst*
dtype0*
value	B :
Z
model_9/weight_1/Pow_1Powmodel_9/weight_1/Pow_1/xmodel_9/weight_1/Pow_1/y*
T0
_
model_9/weight_1/Cast_1Castmodel_9/weight_1/Pow_1*

DstT0*

SrcT0*
Truncate( 
_
model_9/weight_1/truediv_2RealDivmodel_9/weight_1/mul_4model_9/weight_1/Cast_1*
T0
B
model_9/weight_1/Neg_2Negmodel_9/weight_1/truediv_2*
T0
F
model_9/weight_1/Round_1Roundmodel_9/weight_1/truediv_2*
T0
Z
model_9/weight_1/add_4AddV2model_9/weight_1/Neg_2model_9/weight_1/Round_1*
T0
P
model_9/weight_1/StopGradient_2StopGradientmodel_9/weight_1/add_4*
T0
e
model_9/weight_1/add_5AddV2model_9/weight_1/truediv_2model_9/weight_1/StopGradient_2*
T0
W
*model_9/weight_1/clip_by_value_1/Minimum/yConst*
dtype0*
valueB
 * ��C
�
(model_9/weight_1/clip_by_value_1/MinimumMinimummodel_9/weight_1/add_5*model_9/weight_1/clip_by_value_1/Minimum/y*
T0
O
"model_9/weight_1/clip_by_value_1/yConst*
dtype0*
valueB
 *   �
�
 model_9/weight_1/clip_by_value_1Maximum(model_9/weight_1/clip_by_value_1/Minimum"model_9/weight_1/clip_by_value_1/y*
T0
a
model_9/weight_1/mul_5Mulmodel_9/weight_1/Cast_1 model_9/weight_1/clip_by_value_1*
T0
I
model_9/weight_1/truediv_3/yConst*
dtype0*
valueB
 *   D
d
model_9/weight_1/truediv_3RealDivmodel_9/weight_1/mul_5model_9/weight_1/truediv_3/y*
T0
\
model_9/weight_1/mul_6Mulmodel_9/weight_1/mul_6/xmodel_9/weight_1/truediv_3*
T0
X
model_9/weight_1/add_6AddV2model_9/weight_1/Neg_3model_9/weight_1/mul_6*
T0
X
model_9/weight_1/mul_7Mulmodel_9/weight_1/mul_7/xmodel_9/weight_1/add_6*
T0
P
model_9/weight_1/StopGradient_3StopGradientmodel_9/weight_1/mul_7*
T0
l
model_9/weight_1/add_7AddV2!model_9/weight_1/ReadVariableOp_5model_9/weight_1/StopGradient_3*
T0
t
model_9/weight_1/BiasAddBiasAddmodel_9/weight_1/MatMulmodel_9/weight_1/add_7*
T0*
data_formatNHWC
H
model_9/q_activation_6/Pow_1/xConst*
dtype0*
value	B :
H
model_9/q_activation_6/Pow_1/yConst*
dtype0*
value	B :
l
model_9/q_activation_6/Pow_1Powmodel_9/q_activation_6/Pow_1/xmodel_9/q_activation_6/Pow_1/y*
T0
k
model_9/q_activation_6/Cast_1Castmodel_9/q_activation_6/Pow_1*

DstT0*

SrcT0*
Truncate( 
I
model_9/q_activation_6/ConstConst*
dtype0*
valueB
 *   @
I
model_9/q_activation_6/Cast_2/xConst*
dtype0*
value	B :
n
model_9/q_activation_6/Cast_2Castmodel_9/q_activation_6/Cast_2/x*

DstT0*

SrcT0*
Truncate( 
I
model_9/q_activation_6/sub/yConst*
dtype0*
valueB
 *  �A
g
model_9/q_activation_6/subSubmodel_9/q_activation_6/Cast_2model_9/q_activation_6/sub/y*
T0
f
model_9/q_activation_6/Pow_2Powmodel_9/q_activation_6/Constmodel_9/q_activation_6/sub*
T0
i
model_9/q_activation_6/sub_1Submodel_9/q_activation_6/Cast_1model_9/q_activation_6/Pow_2*
T0
n
 model_9/q_activation_6/LessEqual	LessEqualmodel_9/weight_1/BiasAddmodel_9/q_activation_6/sub_1*
T0
F
model_9/q_activation_6/ReluRelumodel_9/weight_1/BiasAdd*
T0
b
&model_9/q_activation_6/ones_like/ShapeShapemodel_9/weight_1/BiasAdd*
T0*
out_type0
S
&model_9/q_activation_6/ones_like/ConstConst*
dtype0*
valueB
 *  �?
�
 model_9/q_activation_6/ones_likeFill&model_9/q_activation_6/ones_like/Shape&model_9/q_activation_6/ones_like/Const*
T0*

index_type0
i
model_9/q_activation_6/sub_2Submodel_9/q_activation_6/Cast_1model_9/q_activation_6/Pow_2*
T0
j
model_9/q_activation_6/mulMul model_9/q_activation_6/ones_likemodel_9/q_activation_6/sub_2*
T0
�
model_9/q_activation_6/SelectV2SelectV2 model_9/q_activation_6/LessEqualmodel_9/q_activation_6/Relumodel_9/q_activation_6/mul*
T0
M
model_9/q_activation_6/Neg_1Negmodel_9/q_activation_6/SelectV2*
T0
F
model_9/q_activation_6/Pow/xConst*
dtype0*
value	B :
F
model_9/q_activation_6/Pow/yConst*
dtype0*
value	B :
f
model_9/q_activation_6/PowPowmodel_9/q_activation_6/Pow/xmodel_9/q_activation_6/Pow/y*
T0
g
model_9/q_activation_6/CastCastmodel_9/q_activation_6/Pow*

DstT0*

SrcT0*
Truncate( 
c
model_9/q_activation_6/mul_1Mulmodel_9/weight_1/BiasAddmodel_9/q_activation_6/Cast*
T0
o
model_9/q_activation_6/truedivRealDivmodel_9/q_activation_6/mul_1model_9/q_activation_6/Cast_1*
T0
J
model_9/q_activation_6/NegNegmodel_9/q_activation_6/truediv*
T0
N
model_9/q_activation_6/RoundRoundmodel_9/q_activation_6/truediv*
T0
f
model_9/q_activation_6/addAddV2model_9/q_activation_6/Negmodel_9/q_activation_6/Round*
T0
X
#model_9/q_activation_6/StopGradientStopGradientmodel_9/q_activation_6/add*
T0
s
model_9/q_activation_6/add_1AddV2model_9/q_activation_6/truediv#model_9/q_activation_6/StopGradient*
T0
o
 model_9/q_activation_6/truediv_1RealDivmodel_9/q_activation_6/add_1model_9/q_activation_6/Cast*
T0
K
model_9/q_activation_6/sub_3/xConst*
dtype0*
valueB
 *  �?
O
"model_9/q_activation_6/truediv_2/xConst*
dtype0*
valueB
 *  �?
u
 model_9/q_activation_6/truediv_2RealDiv"model_9/q_activation_6/truediv_2/xmodel_9/q_activation_6/Cast*
T0
n
model_9/q_activation_6/sub_3Submodel_9/q_activation_6/sub_3/x model_9/q_activation_6/truediv_2*
T0
�
,model_9/q_activation_6/clip_by_value/MinimumMinimum model_9/q_activation_6/truediv_1model_9/q_activation_6/sub_3*
T0
S
&model_9/q_activation_6/clip_by_value/yConst*
dtype0*
valueB
 *    
�
$model_9/q_activation_6/clip_by_valueMaximum,model_9/q_activation_6/clip_by_value/Minimum&model_9/q_activation_6/clip_by_value/y*
T0
q
model_9/q_activation_6/mul_2Mulmodel_9/q_activation_6/Cast_1$model_9/q_activation_6/clip_by_value*
T0
j
model_9/q_activation_6/add_2AddV2model_9/q_activation_6/Neg_1model_9/q_activation_6/mul_2*
T0
j
model_9/q_activation_6/mul_3Mulmodel_9/q_activation_6/mul_3/xmodel_9/q_activation_6/add_2*
T0
\
%model_9/q_activation_6/StopGradient_1StopGradientmodel_9/q_activation_6/mul_3*
T0
v
model_9/q_activation_6/add_3AddV2model_9/q_activation_6/SelectV2%model_9/q_activation_6/StopGradient_1*
T0
�
(model_9/weight_2/ReadVariableOp/resourceConst*
dtype0*�
value�B�

"�<��=[�>    �i>�M��ȏ�R� �    	�>    �L->$bξX'���S>�>�5G>]ť>   �t)�>    ���>�U�>�����B�>|Y >    �?   �T��    Ss־��L>   �   �   �           �o��               �   �   �   �                �/���Ӿ�\)��>   �o�1>V���    �v��   �           �   �           �   �   �       �   �   �   �           �   �        K�<>C&_�)�e>�7��   �P����=    �襾   ����=L�d> ��>�o?<z�>`�>!��   �G��   �
`
!model_9/weight_2/ReadVariableOp_2Identity(model_9/weight_2/ReadVariableOp/resource*
T0
E
model_9/weight_2/mul_3/xConst*
dtype0*
valueB
 *  �?
`
!model_9/weight_2/ReadVariableOp_1Identity(model_9/weight_2/ReadVariableOp/resource*
T0
I
model_9/weight_2/Neg_1Neg!model_9/weight_2/ReadVariableOp_1*
T0
E
model_9/weight_2/mul_2/xConst*
dtype0*
valueB
 *  �?
^
model_9/weight_2/ReadVariableOpIdentity(model_9/weight_2/ReadVariableOp/resource*
T0
C
model_9/weight_2/mul/yConst*
dtype0*
valueB
 *   D
]
model_9/weight_2/mulMulmodel_9/weight_2/ReadVariableOpmodel_9/weight_2/mul/y*
T0
@
model_9/weight_2/Pow/xConst*
dtype0*
value	B :
@
model_9/weight_2/Pow/yConst*
dtype0*
value	B :
T
model_9/weight_2/PowPowmodel_9/weight_2/Pow/xmodel_9/weight_2/Pow/y*
T0
[
model_9/weight_2/CastCastmodel_9/weight_2/Pow*

DstT0*

SrcT0*
Truncate( 
Y
model_9/weight_2/truedivRealDivmodel_9/weight_2/mulmodel_9/weight_2/Cast*
T0
>
model_9/weight_2/NegNegmodel_9/weight_2/truediv*
T0
B
model_9/weight_2/RoundRoundmodel_9/weight_2/truediv*
T0
T
model_9/weight_2/addAddV2model_9/weight_2/Negmodel_9/weight_2/Round*
T0
L
model_9/weight_2/StopGradientStopGradientmodel_9/weight_2/add*
T0
a
model_9/weight_2/add_1AddV2model_9/weight_2/truedivmodel_9/weight_2/StopGradient*
T0
U
(model_9/weight_2/clip_by_value/Minimum/yConst*
dtype0*
valueB
 * ��C
|
&model_9/weight_2/clip_by_value/MinimumMinimummodel_9/weight_2/add_1(model_9/weight_2/clip_by_value/Minimum/y*
T0
M
 model_9/weight_2/clip_by_value/yConst*
dtype0*
valueB
 *   �
|
model_9/weight_2/clip_by_valueMaximum&model_9/weight_2/clip_by_value/Minimum model_9/weight_2/clip_by_value/y*
T0
]
model_9/weight_2/mul_1Mulmodel_9/weight_2/Castmodel_9/weight_2/clip_by_value*
T0
I
model_9/weight_2/truediv_1/yConst*
dtype0*
valueB
 *   D
d
model_9/weight_2/truediv_1RealDivmodel_9/weight_2/mul_1model_9/weight_2/truediv_1/y*
T0
\
model_9/weight_2/mul_2Mulmodel_9/weight_2/mul_2/xmodel_9/weight_2/truediv_1*
T0
X
model_9/weight_2/add_2AddV2model_9/weight_2/Neg_1model_9/weight_2/mul_2*
T0
X
model_9/weight_2/mul_3Mulmodel_9/weight_2/mul_3/xmodel_9/weight_2/add_2*
T0
P
model_9/weight_2/StopGradient_1StopGradientmodel_9/weight_2/mul_3*
T0
l
model_9/weight_2/add_3AddV2!model_9/weight_2/ReadVariableOp_2model_9/weight_2/StopGradient_1*
T0
�
model_9/weight_2/MatMulMatMulmodel_9/q_activation_6/add_3model_9/weight_2/add_3*
T0*
transpose_a( *
transpose_b( 

*model_9/weight_2/ReadVariableOp_3/resourceConst*
dtype0*=
value4B2
"(x8J��*=�/R�3��m��z�;���    �}��    
b
!model_9/weight_2/ReadVariableOp_5Identity*model_9/weight_2/ReadVariableOp_3/resource*
T0
E
model_9/weight_2/mul_7/xConst*
dtype0*
valueB
 *  �?
b
!model_9/weight_2/ReadVariableOp_4Identity*model_9/weight_2/ReadVariableOp_3/resource*
T0
I
model_9/weight_2/Neg_3Neg!model_9/weight_2/ReadVariableOp_4*
T0
E
model_9/weight_2/mul_6/xConst*
dtype0*
valueB
 *  �?
b
!model_9/weight_2/ReadVariableOp_3Identity*model_9/weight_2/ReadVariableOp_3/resource*
T0
E
model_9/weight_2/mul_4/yConst*
dtype0*
valueB
 *   D
c
model_9/weight_2/mul_4Mul!model_9/weight_2/ReadVariableOp_3model_9/weight_2/mul_4/y*
T0
B
model_9/weight_2/Pow_1/xConst*
dtype0*
value	B :
B
model_9/weight_2/Pow_1/yConst*
dtype0*
value	B :
Z
model_9/weight_2/Pow_1Powmodel_9/weight_2/Pow_1/xmodel_9/weight_2/Pow_1/y*
T0
_
model_9/weight_2/Cast_1Castmodel_9/weight_2/Pow_1*

DstT0*

SrcT0*
Truncate( 
_
model_9/weight_2/truediv_2RealDivmodel_9/weight_2/mul_4model_9/weight_2/Cast_1*
T0
B
model_9/weight_2/Neg_2Negmodel_9/weight_2/truediv_2*
T0
F
model_9/weight_2/Round_1Roundmodel_9/weight_2/truediv_2*
T0
Z
model_9/weight_2/add_4AddV2model_9/weight_2/Neg_2model_9/weight_2/Round_1*
T0
P
model_9/weight_2/StopGradient_2StopGradientmodel_9/weight_2/add_4*
T0
e
model_9/weight_2/add_5AddV2model_9/weight_2/truediv_2model_9/weight_2/StopGradient_2*
T0
W
*model_9/weight_2/clip_by_value_1/Minimum/yConst*
dtype0*
valueB
 * ��C
�
(model_9/weight_2/clip_by_value_1/MinimumMinimummodel_9/weight_2/add_5*model_9/weight_2/clip_by_value_1/Minimum/y*
T0
O
"model_9/weight_2/clip_by_value_1/yConst*
dtype0*
valueB
 *   �
�
 model_9/weight_2/clip_by_value_1Maximum(model_9/weight_2/clip_by_value_1/Minimum"model_9/weight_2/clip_by_value_1/y*
T0
a
model_9/weight_2/mul_5Mulmodel_9/weight_2/Cast_1 model_9/weight_2/clip_by_value_1*
T0
I
model_9/weight_2/truediv_3/yConst*
dtype0*
valueB
 *   D
d
model_9/weight_2/truediv_3RealDivmodel_9/weight_2/mul_5model_9/weight_2/truediv_3/y*
T0
\
model_9/weight_2/mul_6Mulmodel_9/weight_2/mul_6/xmodel_9/weight_2/truediv_3*
T0
X
model_9/weight_2/add_6AddV2model_9/weight_2/Neg_3model_9/weight_2/mul_6*
T0
X
model_9/weight_2/mul_7Mulmodel_9/weight_2/mul_7/xmodel_9/weight_2/add_6*
T0
P
model_9/weight_2/StopGradient_3StopGradientmodel_9/weight_2/mul_7*
T0
l
model_9/weight_2/add_7AddV2!model_9/weight_2/ReadVariableOp_5model_9/weight_2/StopGradient_3*
T0
t
model_9/weight_2/BiasAddBiasAddmodel_9/weight_2/MatMulmodel_9/weight_2/add_7*
T0*
data_formatNHWC
H
model_9/q_activation_7/Pow_1/xConst*
dtype0*
value	B :
H
model_9/q_activation_7/Pow_1/yConst*
dtype0*
value	B :
l
model_9/q_activation_7/Pow_1Powmodel_9/q_activation_7/Pow_1/xmodel_9/q_activation_7/Pow_1/y*
T0
k
model_9/q_activation_7/Cast_1Castmodel_9/q_activation_7/Pow_1*

DstT0*

SrcT0*
Truncate( 
I
model_9/q_activation_7/ConstConst*
dtype0*
valueB
 *   @
I
model_9/q_activation_7/Cast_2/xConst*
dtype0*
value	B :
n
model_9/q_activation_7/Cast_2Castmodel_9/q_activation_7/Cast_2/x*

DstT0*

SrcT0*
Truncate( 
I
model_9/q_activation_7/sub/yConst*
dtype0*
valueB
 *  �A
g
model_9/q_activation_7/subSubmodel_9/q_activation_7/Cast_2model_9/q_activation_7/sub/y*
T0
f
model_9/q_activation_7/Pow_2Powmodel_9/q_activation_7/Constmodel_9/q_activation_7/sub*
T0
i
model_9/q_activation_7/sub_1Submodel_9/q_activation_7/Cast_1model_9/q_activation_7/Pow_2*
T0
n
 model_9/q_activation_7/LessEqual	LessEqualmodel_9/weight_2/BiasAddmodel_9/q_activation_7/sub_1*
T0
F
model_9/q_activation_7/ReluRelumodel_9/weight_2/BiasAdd*
T0
b
&model_9/q_activation_7/ones_like/ShapeShapemodel_9/weight_2/BiasAdd*
T0*
out_type0
S
&model_9/q_activation_7/ones_like/ConstConst*
dtype0*
valueB
 *  �?
�
 model_9/q_activation_7/ones_likeFill&model_9/q_activation_7/ones_like/Shape&model_9/q_activation_7/ones_like/Const*
T0*

index_type0
i
model_9/q_activation_7/sub_2Submodel_9/q_activation_7/Cast_1model_9/q_activation_7/Pow_2*
T0
j
model_9/q_activation_7/mulMul model_9/q_activation_7/ones_likemodel_9/q_activation_7/sub_2*
T0
�
model_9/q_activation_7/SelectV2SelectV2 model_9/q_activation_7/LessEqualmodel_9/q_activation_7/Relumodel_9/q_activation_7/mul*
T0
M
model_9/q_activation_7/Neg_1Negmodel_9/q_activation_7/SelectV2*
T0
F
model_9/q_activation_7/Pow/xConst*
dtype0*
value	B :
F
model_9/q_activation_7/Pow/yConst*
dtype0*
value	B :
f
model_9/q_activation_7/PowPowmodel_9/q_activation_7/Pow/xmodel_9/q_activation_7/Pow/y*
T0
g
model_9/q_activation_7/CastCastmodel_9/q_activation_7/Pow*

DstT0*

SrcT0*
Truncate( 
c
model_9/q_activation_7/mul_1Mulmodel_9/weight_2/BiasAddmodel_9/q_activation_7/Cast*
T0
o
model_9/q_activation_7/truedivRealDivmodel_9/q_activation_7/mul_1model_9/q_activation_7/Cast_1*
T0
J
model_9/q_activation_7/NegNegmodel_9/q_activation_7/truediv*
T0
N
model_9/q_activation_7/RoundRoundmodel_9/q_activation_7/truediv*
T0
f
model_9/q_activation_7/addAddV2model_9/q_activation_7/Negmodel_9/q_activation_7/Round*
T0
X
#model_9/q_activation_7/StopGradientStopGradientmodel_9/q_activation_7/add*
T0
s
model_9/q_activation_7/add_1AddV2model_9/q_activation_7/truediv#model_9/q_activation_7/StopGradient*
T0
o
 model_9/q_activation_7/truediv_1RealDivmodel_9/q_activation_7/add_1model_9/q_activation_7/Cast*
T0
K
model_9/q_activation_7/sub_3/xConst*
dtype0*
valueB
 *  �?
O
"model_9/q_activation_7/truediv_2/xConst*
dtype0*
valueB
 *  �?
u
 model_9/q_activation_7/truediv_2RealDiv"model_9/q_activation_7/truediv_2/xmodel_9/q_activation_7/Cast*
T0
n
model_9/q_activation_7/sub_3Submodel_9/q_activation_7/sub_3/x model_9/q_activation_7/truediv_2*
T0
�
,model_9/q_activation_7/clip_by_value/MinimumMinimum model_9/q_activation_7/truediv_1model_9/q_activation_7/sub_3*
T0
S
&model_9/q_activation_7/clip_by_value/yConst*
dtype0*
valueB
 *    
�
$model_9/q_activation_7/clip_by_valueMaximum,model_9/q_activation_7/clip_by_value/Minimum&model_9/q_activation_7/clip_by_value/y*
T0
q
model_9/q_activation_7/mul_2Mulmodel_9/q_activation_7/Cast_1$model_9/q_activation_7/clip_by_value*
T0
j
model_9/q_activation_7/add_2AddV2model_9/q_activation_7/Neg_1model_9/q_activation_7/mul_2*
T0
j
model_9/q_activation_7/mul_3Mulmodel_9/q_activation_7/mul_3/xmodel_9/q_activation_7/add_2*
T0
\
%model_9/q_activation_7/StopGradient_1StopGradientmodel_9/q_activation_7/mul_3*
T0
v
model_9/q_activation_7/add_3AddV2model_9/q_activation_7/SelectV2%model_9/q_activation_7/StopGradient_1*
T0
�
,model_9/weight_final/ReadVariableOp/resourceConst*
dtype0*A
value8B6
"(�߻��S�}�>��ý    E7�>.�-=   ����>    
h
%model_9/weight_final/ReadVariableOp_2Identity,model_9/weight_final/ReadVariableOp/resource*
T0
I
model_9/weight_final/mul_3/xConst*
dtype0*
valueB
 *  �?
h
%model_9/weight_final/ReadVariableOp_1Identity,model_9/weight_final/ReadVariableOp/resource*
T0
Q
model_9/weight_final/Neg_1Neg%model_9/weight_final/ReadVariableOp_1*
T0
I
model_9/weight_final/mul_2/xConst*
dtype0*
valueB
 *  �?
f
#model_9/weight_final/ReadVariableOpIdentity,model_9/weight_final/ReadVariableOp/resource*
T0
G
model_9/weight_final/mul/yConst*
dtype0*
valueB
 *   D
i
model_9/weight_final/mulMul#model_9/weight_final/ReadVariableOpmodel_9/weight_final/mul/y*
T0
D
model_9/weight_final/Pow/xConst*
dtype0*
value	B :
D
model_9/weight_final/Pow/yConst*
dtype0*
value	B :
`
model_9/weight_final/PowPowmodel_9/weight_final/Pow/xmodel_9/weight_final/Pow/y*
T0
c
model_9/weight_final/CastCastmodel_9/weight_final/Pow*

DstT0*

SrcT0*
Truncate( 
e
model_9/weight_final/truedivRealDivmodel_9/weight_final/mulmodel_9/weight_final/Cast*
T0
F
model_9/weight_final/NegNegmodel_9/weight_final/truediv*
T0
J
model_9/weight_final/RoundRoundmodel_9/weight_final/truediv*
T0
`
model_9/weight_final/addAddV2model_9/weight_final/Negmodel_9/weight_final/Round*
T0
T
!model_9/weight_final/StopGradientStopGradientmodel_9/weight_final/add*
T0
m
model_9/weight_final/add_1AddV2model_9/weight_final/truediv!model_9/weight_final/StopGradient*
T0
Y
,model_9/weight_final/clip_by_value/Minimum/yConst*
dtype0*
valueB
 * ��C
�
*model_9/weight_final/clip_by_value/MinimumMinimummodel_9/weight_final/add_1,model_9/weight_final/clip_by_value/Minimum/y*
T0
Q
$model_9/weight_final/clip_by_value/yConst*
dtype0*
valueB
 *   �
�
"model_9/weight_final/clip_by_valueMaximum*model_9/weight_final/clip_by_value/Minimum$model_9/weight_final/clip_by_value/y*
T0
i
model_9/weight_final/mul_1Mulmodel_9/weight_final/Cast"model_9/weight_final/clip_by_value*
T0
M
 model_9/weight_final/truediv_1/yConst*
dtype0*
valueB
 *   D
p
model_9/weight_final/truediv_1RealDivmodel_9/weight_final/mul_1 model_9/weight_final/truediv_1/y*
T0
h
model_9/weight_final/mul_2Mulmodel_9/weight_final/mul_2/xmodel_9/weight_final/truediv_1*
T0
d
model_9/weight_final/add_2AddV2model_9/weight_final/Neg_1model_9/weight_final/mul_2*
T0
d
model_9/weight_final/mul_3Mulmodel_9/weight_final/mul_3/xmodel_9/weight_final/add_2*
T0
X
#model_9/weight_final/StopGradient_1StopGradientmodel_9/weight_final/mul_3*
T0
x
model_9/weight_final/add_3AddV2%model_9/weight_final/ReadVariableOp_2#model_9/weight_final/StopGradient_1*
T0
�
model_9/weight_final/MatMulMatMulmodel_9/q_activation_7/add_3model_9/weight_final/add_3*
T0*
transpose_a( *
transpose_b( 
_
.model_9/weight_final/ReadVariableOp_3/resourceConst*
dtype0*
valueB*<�Ļ
j
%model_9/weight_final/ReadVariableOp_5Identity.model_9/weight_final/ReadVariableOp_3/resource*
T0
I
model_9/weight_final/mul_7/xConst*
dtype0*
valueB
 *  �?
j
%model_9/weight_final/ReadVariableOp_4Identity.model_9/weight_final/ReadVariableOp_3/resource*
T0
Q
model_9/weight_final/Neg_3Neg%model_9/weight_final/ReadVariableOp_4*
T0
I
model_9/weight_final/mul_6/xConst*
dtype0*
valueB
 *  �?
j
%model_9/weight_final/ReadVariableOp_3Identity.model_9/weight_final/ReadVariableOp_3/resource*
T0
I
model_9/weight_final/mul_4/yConst*
dtype0*
valueB
 *   D
o
model_9/weight_final/mul_4Mul%model_9/weight_final/ReadVariableOp_3model_9/weight_final/mul_4/y*
T0
F
model_9/weight_final/Pow_1/xConst*
dtype0*
value	B :
F
model_9/weight_final/Pow_1/yConst*
dtype0*
value	B :
f
model_9/weight_final/Pow_1Powmodel_9/weight_final/Pow_1/xmodel_9/weight_final/Pow_1/y*
T0
g
model_9/weight_final/Cast_1Castmodel_9/weight_final/Pow_1*

DstT0*

SrcT0*
Truncate( 
k
model_9/weight_final/truediv_2RealDivmodel_9/weight_final/mul_4model_9/weight_final/Cast_1*
T0
J
model_9/weight_final/Neg_2Negmodel_9/weight_final/truediv_2*
T0
N
model_9/weight_final/Round_1Roundmodel_9/weight_final/truediv_2*
T0
f
model_9/weight_final/add_4AddV2model_9/weight_final/Neg_2model_9/weight_final/Round_1*
T0
X
#model_9/weight_final/StopGradient_2StopGradientmodel_9/weight_final/add_4*
T0
q
model_9/weight_final/add_5AddV2model_9/weight_final/truediv_2#model_9/weight_final/StopGradient_2*
T0
[
.model_9/weight_final/clip_by_value_1/Minimum/yConst*
dtype0*
valueB
 * ��C
�
,model_9/weight_final/clip_by_value_1/MinimumMinimummodel_9/weight_final/add_5.model_9/weight_final/clip_by_value_1/Minimum/y*
T0
S
&model_9/weight_final/clip_by_value_1/yConst*
dtype0*
valueB
 *   �
�
$model_9/weight_final/clip_by_value_1Maximum,model_9/weight_final/clip_by_value_1/Minimum&model_9/weight_final/clip_by_value_1/y*
T0
m
model_9/weight_final/mul_5Mulmodel_9/weight_final/Cast_1$model_9/weight_final/clip_by_value_1*
T0
M
 model_9/weight_final/truediv_3/yConst*
dtype0*
valueB
 *   D
p
model_9/weight_final/truediv_3RealDivmodel_9/weight_final/mul_5 model_9/weight_final/truediv_3/y*
T0
h
model_9/weight_final/mul_6Mulmodel_9/weight_final/mul_6/xmodel_9/weight_final/truediv_3*
T0
d
model_9/weight_final/add_6AddV2model_9/weight_final/Neg_3model_9/weight_final/mul_6*
T0
d
model_9/weight_final/mul_7Mulmodel_9/weight_final/mul_7/xmodel_9/weight_final/add_6*
T0
X
#model_9/weight_final/StopGradient_3StopGradientmodel_9/weight_final/mul_7*
T0
x
model_9/weight_final/add_7AddV2%model_9/weight_final/ReadVariableOp_5#model_9/weight_final/StopGradient_3*
T0
�
model_9/weight_final/BiasAddBiasAddmodel_9/weight_final/MatMulmodel_9/weight_final/add_7*
T0*
data_formatNHWC
H
model_9/q_activation_8/Pow_1/xConst*
dtype0*
value	B :
H
model_9/q_activation_8/Pow_1/yConst*
dtype0*
value	B :
l
model_9/q_activation_8/Pow_1Powmodel_9/q_activation_8/Pow_1/xmodel_9/q_activation_8/Pow_1/y*
T0
k
model_9/q_activation_8/Cast_1Castmodel_9/q_activation_8/Pow_1*

DstT0*

SrcT0*
Truncate( 
I
model_9/q_activation_8/ConstConst*
dtype0*
valueB
 *   @
I
model_9/q_activation_8/Cast_2/xConst*
dtype0*
value	B :
n
model_9/q_activation_8/Cast_2Castmodel_9/q_activation_8/Cast_2/x*

DstT0*

SrcT0*
Truncate( 
I
model_9/q_activation_8/sub/yConst*
dtype0*
valueB
 *  �A
g
model_9/q_activation_8/subSubmodel_9/q_activation_8/Cast_2model_9/q_activation_8/sub/y*
T0
f
model_9/q_activation_8/Pow_2Powmodel_9/q_activation_8/Constmodel_9/q_activation_8/sub*
T0
i
model_9/q_activation_8/sub_1Submodel_9/q_activation_8/Cast_1model_9/q_activation_8/Pow_2*
T0
r
 model_9/q_activation_8/LessEqual	LessEqualmodel_9/weight_final/BiasAddmodel_9/q_activation_8/sub_1*
T0
J
model_9/q_activation_8/ReluRelumodel_9/weight_final/BiasAdd*
T0
f
&model_9/q_activation_8/ones_like/ShapeShapemodel_9/weight_final/BiasAdd*
T0*
out_type0
S
&model_9/q_activation_8/ones_like/ConstConst*
dtype0*
valueB
 *  �?
�
 model_9/q_activation_8/ones_likeFill&model_9/q_activation_8/ones_like/Shape&model_9/q_activation_8/ones_like/Const*
T0*

index_type0
i
model_9/q_activation_8/sub_2Submodel_9/q_activation_8/Cast_1model_9/q_activation_8/Pow_2*
T0
j
model_9/q_activation_8/mulMul model_9/q_activation_8/ones_likemodel_9/q_activation_8/sub_2*
T0
�
model_9/q_activation_8/SelectV2SelectV2 model_9/q_activation_8/LessEqualmodel_9/q_activation_8/Relumodel_9/q_activation_8/mul*
T0
M
model_9/q_activation_8/Neg_1Negmodel_9/q_activation_8/SelectV2*
T0
F
model_9/q_activation_8/Pow/xConst*
dtype0*
value	B :
F
model_9/q_activation_8/Pow/yConst*
dtype0*
value	B :
f
model_9/q_activation_8/PowPowmodel_9/q_activation_8/Pow/xmodel_9/q_activation_8/Pow/y*
T0
g
model_9/q_activation_8/CastCastmodel_9/q_activation_8/Pow*

DstT0*

SrcT0*
Truncate( 
g
model_9/q_activation_8/mul_1Mulmodel_9/weight_final/BiasAddmodel_9/q_activation_8/Cast*
T0
o
model_9/q_activation_8/truedivRealDivmodel_9/q_activation_8/mul_1model_9/q_activation_8/Cast_1*
T0
J
model_9/q_activation_8/NegNegmodel_9/q_activation_8/truediv*
T0
N
model_9/q_activation_8/RoundRoundmodel_9/q_activation_8/truediv*
T0
f
model_9/q_activation_8/addAddV2model_9/q_activation_8/Negmodel_9/q_activation_8/Round*
T0
X
#model_9/q_activation_8/StopGradientStopGradientmodel_9/q_activation_8/add*
T0
s
model_9/q_activation_8/add_1AddV2model_9/q_activation_8/truediv#model_9/q_activation_8/StopGradient*
T0
o
 model_9/q_activation_8/truediv_1RealDivmodel_9/q_activation_8/add_1model_9/q_activation_8/Cast*
T0
K
model_9/q_activation_8/sub_3/xConst*
dtype0*
valueB
 *  �?
O
"model_9/q_activation_8/truediv_2/xConst*
dtype0*
valueB
 *  �?
u
 model_9/q_activation_8/truediv_2RealDiv"model_9/q_activation_8/truediv_2/xmodel_9/q_activation_8/Cast*
T0
n
model_9/q_activation_8/sub_3Submodel_9/q_activation_8/sub_3/x model_9/q_activation_8/truediv_2*
T0
�
,model_9/q_activation_8/clip_by_value/MinimumMinimum model_9/q_activation_8/truediv_1model_9/q_activation_8/sub_3*
T0
S
&model_9/q_activation_8/clip_by_value/yConst*
dtype0*
valueB
 *    
�
$model_9/q_activation_8/clip_by_valueMaximum,model_9/q_activation_8/clip_by_value/Minimum&model_9/q_activation_8/clip_by_value/y*
T0
q
model_9/q_activation_8/mul_2Mulmodel_9/q_activation_8/Cast_1$model_9/q_activation_8/clip_by_value*
T0
j
model_9/q_activation_8/add_2AddV2model_9/q_activation_8/Neg_1model_9/q_activation_8/mul_2*
T0
j
model_9/q_activation_8/mul_3Mulmodel_9/q_activation_8/mul_3/xmodel_9/q_activation_8/add_2*
T0
\
%model_9/q_activation_8/StopGradient_1StopGradientmodel_9/q_activation_8/mul_3*
T0
v
model_9/q_activation_8/add_3AddV2model_9/q_activation_8/SelectV2%model_9/q_activation_8/StopGradient_1*
T0
�
NoOpNoOp ^model_9/weight_1/ReadVariableOp"^model_9/weight_1/ReadVariableOp_1"^model_9/weight_1/ReadVariableOp_2"^model_9/weight_1/ReadVariableOp_3"^model_9/weight_1/ReadVariableOp_4"^model_9/weight_1/ReadVariableOp_5 ^model_9/weight_2/ReadVariableOp"^model_9/weight_2/ReadVariableOp_1"^model_9/weight_2/ReadVariableOp_2"^model_9/weight_2/ReadVariableOp_3"^model_9/weight_2/ReadVariableOp_4"^model_9/weight_2/ReadVariableOp_5$^model_9/weight_final/ReadVariableOp&^model_9/weight_final/ReadVariableOp_1&^model_9/weight_final/ReadVariableOp_2&^model_9/weight_final/ReadVariableOp_3&^model_9/weight_final/ReadVariableOp_4&^model_9/weight_final/ReadVariableOp_5*"
_acd_function_control_output(
B
IdentityIdentitymodel_9/q_activation_8/add_3^NoOp*
T0"�