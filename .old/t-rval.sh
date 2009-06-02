CXXFLAGS+="-DOARG='trace<V>(rosenberg<V>())' -DxOARG='rosenberg<V>()' -D_N=2" make -B  t-rval-ref && t-rval-ref    &&
CXXFLAGS+="-DxOARG='trace<V>(rosenberg<V>())' -DOARG='rosenberg<V>()' -D_N=2" make -B  t-rval-ref && t-rval-ref
