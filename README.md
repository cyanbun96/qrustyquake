# QrustyQuake-Puppy

Good ol' WinQuake! Except SDL3, with some extra options!

---

This fork mainly exists to create a branch for submitting improvments as pull requests, but i'll probably end up with another branch or two with some mangled code to look at in the near future!

# Added Features!

- Fixed crash where particles seen through viewport with 'sizedown' caused lockup.

- Fixed problem of BGM track not looping/repeating during map play.

- Fixed edge case crash where D_SpriteScanEdge could let sprite_spans write past its memory limit.

- Fixed edge case crash where game would randomly freeze due to drawspan overflow. (could not reproduce error, was never constant during testing)

# Planned

- Exploration into overhaul of progs commands, maybe customizing available functions?

- Keeping the organizational structure and readability of the source code top-notch, while not obfuscating much/any of the original code.

- An additional set of custom tools to help provide a workflow for creation of custom games and tech demos, with source code available.

Contributions of any kind are very welcome.

Issue reports are also greatly appreciated.

# Compiling

Linux:

    `cd src && make`

*nix:

    `cd src && make -f Makefile.freebsd/openbsd/...`

Windows: (Compile using Linux)

    `cd src && make -f Makefile.w64`

Native Build:

    `cmake . && make`

# Credits

All credit goes to the original authors of QrustyQuake, their credits list is fully available on their repo page.

Copyright(C) 1996-2001 Id Software, Inc.
Copyright(C) 2002-2009 John Fitzgibbons and others
Copyright(C) 2010-2014 QuakeSpasm developers
Copyright(C) 2024-2025 QrustyQuake developers
