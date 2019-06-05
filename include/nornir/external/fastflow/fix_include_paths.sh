#!/bin/bash
find . -type f -exec sed -i 's#<ff/#<nornir/external/fastflow/ff/#g' {} +