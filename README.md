### Build

```bash
cmake -DCMAKE_BUILD_TYPE=Release -B build
make -j -C build
```

### Launch

```bash
$ cd problems/mixturn
$ PYTHONPATH=../../build/ python3 mixturn.py 
```

### View results

In ParaView open `vtks/scene...vtm`
