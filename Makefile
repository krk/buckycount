# Copyright 2026 Kerem KAT
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


.DEFAULT_GOAL := all
CXX ?= c++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++20 -Wall -Wextra -Wpedantic
PKG_CONFIG ?= pkg-config
PREFIX ?= /usr/local
INSTALL ?= install
TARGET ?= buckycount
SOURCE := .
UNITS := $(addprefix $(SOURCE)/,buckycount.cpp cli.cpp graph_cli.cpp decode.cpp fullerene_rank_total.cpp)

.PHONY: all install clean
all: $(TARGET)
$(TARGET): $(UNITS) $(wildcard $(SOURCE)/*.hpp) Makefile $(EXTRA_DEPS)
	$(PKG_CONFIG) --print-errors --exists flint gmpxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -pthread $(patsubst -I%,-isystem %,$(shell $(PKG_CONFIG) --cflags flint gmpxx)) $(UNITS) $(LDFLAGS) $(shell $(PKG_CONFIG) --libs-only-L flint gmpxx) -lflint -lgmpxx -lgmp -o $@
install: $(TARGET)
	$(INSTALL) -d "$(DESTDIR)$(PREFIX)/bin"
	$(INSTALL) -m 755 "$(TARGET)" "$(DESTDIR)$(PREFIX)/bin/buckycount"
clean:
	$(RM) $(TARGET)
