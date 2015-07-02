#!/bin/env bash

mkdir -p $PREFIX/bin
mkdir -p $PREFIX/lib/QueryFuse

cat <<EOF >$PREFIX/bin/QF_wrapper
#!/bin/env bash
SCRIPT_DIR=\$( cd "\$( dirname "\${BASH_SOURCE[0]}" )" && pwd )
QF_DIR=\$( cd \$SCRIPT_DIR/../lib/QueryFuse && pwd )
python \$QF_DIR/QF_multi_query_wrapper.py \$@
EOF
chmod +x $PREFIX/bin/QF_wrapper
cp ./QueryFuse_v1/* $PREFIX/lib/QueryFuse/
