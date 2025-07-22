# Error in Heasoft Ver.6.35.*
-- Note: in a case of my Mac (Apple M3 Sequoia 15.5) --

There is a bug in a file
```
heasoft-6.35.*/maxi/tasks/mxscancur/mxscancurfunc.c
```
Move to the directory `heasoft-6.35.*/maxi/tasks/mxscancur/`
and replace it with the file [`mxscancurfunc.c`](mxscancurfunc.c) that has been obtained from Heasoft Ver.6.34.

Then, type
```
% hmake clean
% hmake install
```