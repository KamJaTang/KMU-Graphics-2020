#ifndef KMUCS_GRAPHICS_TRANSFORM_HPP
#define KMUCS_GRAPHICS_TRANSFORM_HPP

#include <cmath>
#include "vec.hpp"
#include "mat.hpp"
#include "operator.hpp"

namespace kmuvcl
{
    namespace math
    {
#ifndef M_PI
        const float M_PI = 3.14159265358979323846f;
#endif

        template <typename T>
        mat<4, 4, T> translate(T dx, T dy, T dz)
        {
            mat<4, 4, T> translateMat;

            // TODO: Fill up this function properly

            translateMat[0] = 1;
            translateMat[1] = 0;
            translateMat[2] = 0;
            translateMat[3] = 0;

            translateMat[4] = 0;
            translateMat[5] = 1;
            translateMat[6] = 0;
            translateMat[7] = 0;

            translateMat[8] = 0;
            translateMat[9] = 0;
            translateMat[10] = 1;
            translateMat[11] = 0;

            translateMat[12] = dx;
            translateMat[13] = dy;
            translateMat[14] = dz;
            translateMat[15] = 1;

            return translateMat;
        }

        template <typename T>
        mat<4, 4, T> rotate(T angle, T x, T y, T z)
        {
            mat<4, 4, T> rotateMat;

            // TODO: Fill up this function properly

            float r = angle * M_PI / 180;
            float ux, uy, uz, size;
            size = sqrt(x*x + y*y + z*z);
            ux = x/size;
            uy = y/size;
            uz = z/size;

            rotateMat[0] = cos(r) + ux*ux*(1-cos(r));
            rotateMat[1] = uy*ux*(1-cos(r))+uz*sin(r);
            rotateMat[2] = uz*ux*(1-cos(r))-uy*sin(r);
            rotateMat[3] = 0;

            rotateMat[4] = ux*uy*(1-cos(r))-uz*sin(r);
            rotateMat[5] = cos(r)+uy*uy*(1-cos(r));
            rotateMat[6] = uz*uy*(1-cos(r))+ ux*sin(r);
            rotateMat[7] = 0;

            rotateMat[8] = ux*uz*(1-cos(r))+uy*sin(r);
            rotateMat[9] = uy*uz*(1-cos(r))-ux*sin(r);
            rotateMat[10] = cos(r) + uz*uz*(1-cos(r));
            rotateMat[11] = 0;

            rotateMat[12] = 0;
            rotateMat[13] = 0;
            rotateMat[14] = 0;
            rotateMat[15] = 1;


            return rotateMat;
        }

        template<typename T>
        mat<4, 4, T> scale(T sx, T sy, T sz)
        {
            mat<4, 4, T> scaleMat;

            // TODO: Fill up this function properly

            scaleMat[0] = sx;
            scaleMat[1] = 0;
            scaleMat[2] = 0;
            scaleMat[3] = 0;

            scaleMat[4] = 0;
            scaleMat[5] = sy;
            scaleMat[6] = 0;
            scaleMat[7] = 0;

            scaleMat[8] = 0;
            scaleMat[9] = 0;
            scaleMat[10] = sz;
            scaleMat[11] = 0;

            scaleMat[12] = 0;
            scaleMat[13] = 0;
            scaleMat[14] = 0;
            scaleMat[15] = 1;

            return scaleMat;
        }

        template<typename T>
        mat<4, 4, T> lookAt(T eyeX, T eyeY, T eyeZ, T centerX, T centerY, T centerZ, T upX, T upY, T upZ)
        {
            mat<4, 4, T> viewMat;
            mat<4, 4, T> rMat;
            mat<4, 4, T> tMat;

            vec<3, T> camX;
            vec<3, T> camY;
            vec<3, T> camZ;

            vec<3, T> up;
            up[0] = upX;
            up[1] = upY;
            up[2] = upZ;

            camZ[0] = -centerX + eyeX;
            camZ[1] = -centerY + eyeY;
            camZ[2] = -centerZ + eyeZ;
            float sz = sqrt(camZ[0]*camZ[0] + camZ[1]*camZ[1] + camZ[2]*camZ[2]);
            camZ[0] /=sz;
            camZ[1] /=sz;
            camZ[2] /=sz;

            camX = cross(up, camZ);
            float sx = sqrt(camX[0]*camX[0] + camX[1]*camX[1] + camX[2]*camX[2]);
            camX[0] /= sx;
            camX[1] /= sx;
            camX[2] /= sx;

            camY = cross(camZ, camX);
            float sy = sqrt(camY[0]*camY[0] + camY[1]*camY[1] + camY[2]*camY[2]);
            camY[0] /= sy;
            camY[1] /= sy;
            camY[2] /= sy;

            rMat[0] = camX[0];
            rMat[1] = camY[0];
            rMat[2] = camZ[0];
            rMat[3] = 0;

            rMat[4] = camX[1];
            rMat[5] = camY[1];
            rMat[6] = camZ[1];
            rMat[7] = 0;

            rMat[8] = camX[2];
            rMat[9] = camY[2];
            rMat[10] = camZ[2];
            rMat[11] = 0;

            rMat[12] = 0;
            rMat[13] = 0;
            rMat[14] = 0;
            rMat[15] = 1;


            tMat[0] = 1;
            tMat[1] = 0;
            tMat[2] = 0;
            tMat[3] = 0;

            tMat[4] = 0;
            tMat[5] = 1;
            tMat[6] = 0;
            tMat[7] = 0;

            tMat[8] = 0;
            tMat[9] = 0;
            tMat[10] = 1;
            tMat[11] = 0;

            tMat[12] = -eyeX;
            tMat[13] = -eyeY;
            tMat[14] = -eyeZ;
            tMat[15] = 1;

            viewMat = rMat * tMat;\

            return viewMat;
        }

        template<typename T>
        mat<4, 4, T> ortho(T left, T right, T bottom, T top, T nearVal, T farVal)
        {
            mat<4, 4, T> orthoMat;

            // TODO: Fill up this function properly

            orthoMat[0] = 2/(right-left);
            orthoMat[1] = 0;
            orthoMat[2] = 0;
            orthoMat[3] = 0;

            orthoMat[4] = 0;
            orthoMat[5] = 2/(top-bottom);
            orthoMat[6] = 0;
            orthoMat[7] = 0;

            orthoMat[8] = 0;
            orthoMat[9] = 0;
            orthoMat[10] = -2/(farVal-nearVal);
            orthoMat[11] = 0;

            orthoMat[12] = -(right+left)/(right-left);
            orthoMat[13] = -(top+bottom)/(top-bottom);
            orthoMat[14] = -(farVal+nearVal)/(farVal-nearVal);
            orthoMat[15] = 1;

            return orthoMat;
        }

        template<typename T>
        mat<4, 4, T> frustum(T left, T right, T bottom, T top, T nearVal, T farVal)
        {
           mat<4, 4, T> frustumMat;

           // TODO: Fill up this function properly

           frustumMat[0] = 2*nearVal/(right-left);
           frustumMat[1] = 0;
           frustumMat[2] = 0;
           frustumMat[3] = 0;

           frustumMat[4] = 0;
           frustumMat[5] = 2*nearVal/(top-bottom);
           frustumMat[6] = 0;
           frustumMat[7] = 0;

           frustumMat[8] = (right+left)/(right-left);
           frustumMat[9] = (top+bottom)/(top-bottom);
           frustumMat[10] = -(farVal+nearVal)/(farVal-nearVal);
           frustumMat[11] = -1;

           frustumMat[12] = 0;
           frustumMat[13] = 0;
           frustumMat[14] = -2*farVal*nearVal/(farVal-nearVal);
           frustumMat[15] = 0;

           return frustumMat;
        }

        template<typename T>
        mat<4, 4, T> perspective(T fovy, T aspect, T zNear, T zFar)
        {
          T  right = 0;
          T  top = 0;

          // TODO: Fill up this function properly

          float r = fovy * M_PI / 180 * 1/2;

          top = tan(r) * zNear;
          right = aspect * top;


          return frustum(-right, right, -top, top, zNear, zFar);
        }
    }
}
#endif
