#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define R(mul, shift, x, y)          \
  _ = x;                             \
  x -= mul * y >> shift;             \
  y += mul * _ >> shift;             \
  _ = 3145728 - x * x - y * y >> 11; \
  x = x * _ >> 10;                   \
  y = y * _ >> 10;

int8_t b[1760], z[1760];

void main()
{
  int sA = 1024, cA = 0, sB = 1024, cB = 0, _;
  for (;;)
  {
    memset(b, 32, 1760);  // text buffer
    memset(z, 127, 1760); // z buffer
    int sj = 0, cj = 1024;
    for (int j = 0; j < 90; j++)
    {
      int si = 0, ci = 1024; // sine and cosine of angle i
      for (int i = 0; i < 324; i++)
      {
        int R1 = 1, R2 = 2048, K2 = 5120 * 1024;

        int x0 = R1 * cj + R2,
            x1 = ci * x0 >> 10,
            x2 = cA * sj >> 10,
            x3 = si * x0 >> 10,
            x4 = R1 * x2 - (sA * x3 >> 10),
            x5 = sA * sj >> 10,
            x6 = K2 + R1 * 1024 * x5 + cA * x3,
            x7 = cj * si >> 10,
            x = 40 + 30 * (cB * x1 - sB * x4) / x6,
            y = 12 + 15 * (cB * x4 + sB * x1) / x6,
            N = (-cA * x7 - cB * ((-sA * x7 >> 10) + x2) - ci * (cj * sB >> 10) >> 10) - x5 >> 7;

        int o = x + 80 * y;
        int8_t zz = (x6 - K2) >> 15;
        if (22 > y && y > 0 && x > 0 && 80 > x && zz < z[o])
        {
          z[o] = zz;
          b[o] = ".,-~:;=!*#$@"[N > 0 ? N : 0];
        }
        R(5, 8, ci, si) // rotate i'
      }
      R(9, 7, cj, sj) // rotate j
    }
    for (int k = 0; 1761 > k; k++)
      putchar(k % 80 ? b[k] : 10);
    R(5, 7, cA, sA);
    R(5, 8, cB, sB);
    usleep(15000);
    printf("\x1b[23A");
  }
}
// const float theta_spacing = 0.07;
// const float phi_spacing   = 0.02;

// const float R1 = 1;
// const float R2 = 2;
// const float K2 = 5;
// // Calculate K1 based on screen size: the maximum x-distance occurs
// // roughly at the edge of the torus, which is at x=R1+R2, z=0.  we
// // want that to be displaced 3/8ths of the width of the screen, which
// // is 3/4th of the way from the center to the side of the screen.
// // screen_width*3/8 = K1*(R1+R2)/(K2+0)
// // screen_width*K2*3/(8*(R1+R2)) = K1
// const float K1 = screen_width*K2*3/(8*(R1+R2));

// render_frame(float A, float B) {
//   // precompute sines and cosines of A and B
//   float cosA = cos(A), sinA = sin(A);
//   float cosB = cos(B), sinB = sin(B);

//   char output[0..screen_width, 0..screen_height] = ' ';
//   float zbuffer[0..screen_width, 0..screen_height] = 0;

//   // theta goes around the cross-sectional circle of a torus
//   for (float theta=0; theta < 2*pi; theta += theta_spacing) {
//     // precompute sines and cosines of theta
//     float costheta = cos(theta), sintheta = sin(theta);

//     // phi goes around the center of revolution of a torus
//     for(float phi=0; phi < 2*pi; phi += phi_spacing) {
//       // precompute sines and cosines of phi
//       float cosphi = cos(phi), sinphi = sin(phi);

//       // the x,y coordinate of the circle, before revolving (factored
//       // out of the above equations)
//       float circlex = R2 + R1*costheta;
//       float circley = R1*sintheta;

//       // final 3D (x,y,z) coordinate after rotations, directly from
//       // our math above
//       float x = circlex*(cosB*cosphi + sinA*sinB*sinphi)
//         - circley*cosA*sinB;
//       float y = circlex*(sinB*cosphi - sinA*cosB*sinphi)
//         + circley*cosA*cosB;
//       float z = K2 + cosA*circlex*sinphi + circley*sinA;
//       float ooz = 1/z;  // "one over z"

//       // x and y projection.  note that y is negated here, because y
//       // goes up in 3D space but down on 2D displays.
//       int xp = (int) (screen_width/2 + K1*ooz*x);
//       int yp = (int) (screen_height/2 - K1*ooz*y);

//       // calculate luminance.  ugly, but correct.
//       float L = cosphi*costheta*sinB - cosA*costheta*sinphi -
//         sinA*sintheta + cosB*(cosA*sintheta - costheta*sinA*sinphi);
//       // L ranges from -sqrt(2) to +sqrt(2).  If it's < 0, the surface
//       // is pointing away from us, so we won't bother trying to plot it.
//       if (L > 0) {
//         // test against the z-buffer.  larger 1/z means the pixel is
//         // closer to the viewer than what's already plotted.
//         if(ooz > zbuffer[xp,yp]) {
//           zbuffer[xp, yp] = ooz;
//           int luminance_index = L*8;
//           // luminance_index is now in the range 0..11 (8*sqrt(2) = 11.3)
//           // now we lookup the character corresponding to the
//           // luminance and plot it in our output:
//           output[xp, yp] = ".,-~:;=!*#$@"[luminance_index];
//         }
//       }
//     }
//   }

//   // now, dump output[] to the screen.
//   // bring cursor to "home" location, in just about any currently-used
//   // terminal emulation mode
//   printf("\x1b[H");
//   for (int j = 0; j < screen_height; j++) {
//     for (int i = 0; i < screen_width; i++) {
//       putchar(output[i,j]);
//     }
//     putchar('\n');
//   }

// }
