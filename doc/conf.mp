prologues := 3;
ahangle := 15;
ahlength := 5;
outputtemplate := "%j-%c.eps";

beginfig(0);

a = 35;

b = 5;
z1 = (0,0);
z0 = z1 - (3a, 0);
z2 = z1 + (3a, 0);
z3 = z1;
z4 = z3 + (0, 3.5a);

drawarrow z0 -- z2 withpen pencircle scaled 0.2;
drawarrow z3 -- z4 withpen pencircle scaled 0.2;

label.lrt(btex $x$ etex, z2);
label.rt(btex $t$ etex, z4);

draw (z1 -- z2) rotated 45;
draw (z1 -- z2) rotated 75;
draw (z1 -- z2) rotated 130;

z5 = .9[z1,z2];

draw (z1 -- z5) rotated 132.5;
draw (z1 -- z5) rotated 135;
draw (z1 -- z5) rotated 137.5;

draw (z1 -- z2) rotated 140;

draw (z1 -- z2) rotated 90;

label.rt(btex $SW$ etex, z2 rotated 45);
label.rt(btex $CS$ etex, z2 rotated 75);
label.ulft(btex $RW$ etex, z5 rotated 135);
label.ulft(btex $PJ$ etex, z5 rotated 90);

endfig;

beginfig(1);

z1 = (0,0);
z0 = z1 - (3a, 0);
z2 = z1 + (3a, 0);
z3 = z1;
z4 = z3 + (0, 3.5a);

drawarrow z0 -- z2 withpen pencircle scaled 0.2;
drawarrow z3 -- z4 withpen pencircle scaled 0.2;

label.bot(btex $\xi$ etex, z2);
label.rt(btex $p$ etex, z4);

draw (z1 -- .9z4) shifted (2.5a,0) dashed evenly withpen pencircle scaled 0.2;
draw (z1 -- .9z4) shifted (a,0) dashed evenly withpen pencircle scaled 0.2;

draw (z1 -- .9z4) shifted (-2.2a,0) dashed evenly withpen pencircle scaled 0.2;
draw (z1 -- .9z4) shifted (-2.8a,0) dashed evenly withpen pencircle scaled 0.2;

label.bot(btex $SW$ etex, z1 shifted (2.5a,0));
label.bot(btex $RW$ etex, z1 shifted (-2.5a,0));
label.bot(btex $PJ$ etex, z1);
label.bot(btex $CS$ etex, z1 shifted (a, 0));

draw ((-3a, 0) -- (-2.8a, 0)) shifted .5z4;
draw ((-2.8a, 0){dir -70} .. (-2.2a, -.8a){dir -30}) shifted .5z4;
draw ((-2.2a, -.8a) -- (0, -.8a)) shifted .5z4;
draw ((0, -.8a) -- (2.5a, -.8a)) shifted .9z4;
draw ((2.5a, 0) -- (2.8a, 0)) shifted .5z4;

drawdblarrow (-0.1a, -.8a + .5y4) -- (-0.1a, -.8a + .9y4);
label.lft(btex $\Delta p$ etex, (-0.1a, -.8a + .7y4));

endfig;

beginfig(2);

z1 = (0,0);
z0 = z1 - (3a, 0);
z2 = z1 + (3a, 0);
z3 = z1;
z4 = z3 + (0, 3.5a);

drawarrow z0 -- z2 withpen pencircle scaled 0.2;
drawarrow z3 -- z4 withpen pencircle scaled 0.2;

label.bot(btex $\xi$ etex, z2);
label.rt(btex $\rho$ etex, z4);

draw (z1 -- .9z4) shifted (2.5a,0) dashed evenly withpen pencircle scaled 0.2;
draw (z1 -- .9z4) shifted (a,0) dashed evenly withpen pencircle scaled 0.2;

draw (z1 -- .9z4) shifted (-2.2a,0) dashed evenly withpen pencircle scaled 0.2;
draw (z1 -- .9z4) shifted (-2.8a,0) dashed evenly withpen pencircle scaled 0.2;

label.bot(btex $SW$ etex, z1 shifted (2.5a,0));
label.bot(btex $RW$ etex, z1 shifted (-2.5a,0));
label.bot(btex $PJ$ etex, z1);
label.bot(btex $CS$ etex, z1 shifted (a, 0));

draw ((-3a, 0) -- (-2.8a, 0)) shifted .5z4;
draw ((-2.8a, 0){dir -70} .. (-2.2a, -.5a){dir -10}) shifted .5z4;
draw ((-2.2a, -.5a) -- (a, -.5a)) shifted .5z4;
draw ((a, -.8a) -- (2.5a, -.8a)) shifted .8z4;
draw ((2.5a, 0) -- (2.8a, 0)) shifted .5z4;

endfig;

beginfig(3);

z1 = (0,0);
z0 = z1 - (3a, 0);
z2 = z1 + (3a, 0);
z3 = z1;
z4 = z3 + (0, 3.5a);

drawarrow z0 -- z2 withpen pencircle scaled 0.2;
drawarrow z3 -- z4 withpen pencircle scaled 0.2;

label.bot(btex $\xi$ etex, z2);
label.rt(btex $u$ etex, z4);

draw (z1 -- .9z4) shifted (2.5a,0) dashed evenly withpen pencircle scaled 0.2;
draw (z1 -- .9z4) shifted (a,0) dashed evenly withpen pencircle scaled 0.2;

draw (z1 -- .9z4) shifted (-2.2a,0) dashed evenly withpen pencircle scaled 0.2;
draw (z1 -- .9z4) shifted (-2.8a,0) dashed evenly withpen pencircle scaled 0.2;

label.bot(btex $SW$ etex, z1 shifted (2.5a,0));
label.bot(btex $RW$ etex, z1 shifted (-2.5a,0));
label.bot(btex $PJ$ etex, z1);
label.bot(btex $CS$ etex, z1 shifted (a, 0));

draw ((-3a, 0) -- (-2.8a, 0)) shifted .3z4;
draw ((-2.8a, 0) -- (-2.2a, a)) shifted .3z4;
draw ((-2.2a, a) -- (2.5a, a)) shifted .3z4;
draw ((2.5a, 0) -- (2.8a, 0)) shifted .3z4;

label.ulft(btex $u_0$ etex, (0, a + .3y4));

endfig;

beginfig(4);

z1 = (0,0);
z0 = z1 - (3a, 0);
z2 = z1 + (3a, 0);
z3 = z1;
z4 = z3 + (0, 3.5a);

drawarrow z0 -- z2 withpen pencircle scaled 0.2;
drawarrow z3 -- z4 withpen pencircle scaled 0.2;

label.lrt(btex $u$ etex, z2);
label.rt(btex $p$ etex, z4);

z10 = (-2a, 2.5a);
z11 = (-1.5a, .5a);
z12 = (a, 0);
z13 = (-2.8a, 3a);
z14 = (2.8a, 1.2a);
z15 = (-2.8a, .4a);
z16 = (2.8a, 1.5a);

dotlabel.urt(btex $(u_1, p_1)$ etex, z10);
dotlabel.ulft(btex $(u_2, p_2)$ etex, z11);

draw (z1 -- .9z4) shifted z12 dashed evenly withpen pencircle scaled 0.2;
label.bot(btex $u_0$ etex, z12);

path lef, rig;

lef = z13 ... z10 ... z14;
rig = z15 ... z11 ... z16;

draw lef;
draw rig;

z17 = ((x12, 0) -- (x12, 10a)) intersectionpoint lef;
z18 = ((x12, 0) -- (x12, 10a)) intersectionpoint rig;

dotlabel.urt(btex $(u_0, P_1)$ etex, z17);
dotlabel.lrt(btex $(u_0, P_2)$ etex, z18);

drawdblarrow (z17 -- z18) shifted (-.1a, 0);
label.lft(btex $\Delta p$ etex, .5[z17,z18] shifted (-.1a, 0));

endfig;

end
