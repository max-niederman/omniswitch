{
  description = "omniswitch — keyswitch with electronically simulated force curves (solenoid + magnet + position sensor)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" ];
      # allowUnfree for the FEMM binary distribution (Aladdin Free Public License)
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f (import nixpkgs {
        inherit system;
        config.allowUnfree = true;
      }));
    in
    {
      packages = forAllSystems (pkgs: rec {
        # MSVC 2008 SP1 x64 runtime (mfc90/msvcr90/msvcp90): femm.exe links MFC 9.0
        # dynamically and Wine doesn't provide it. Extracted natively with cabextract
        # and deployed app-local next to femm.exe.
        vc90-runtime = pkgs.stdenvNoCC.mkDerivation {
          pname = "vc90-runtime-x64";
          version = "9.0.30729.6161";
          src = pkgs.fetchurl {
            url = "https://download.microsoft.com/download/5/D/8/5D8C65CB-C849-4025-8E95-C3966CAFD8AE/vcredist_x64.exe";
            hash = "sha256-xeJzpKFqtNVHHpHHR3cZovRd2tt2x/mKOPpQdKaDhlQ=";
          };
          dontUnpack = true;
          nativeBuildInputs = [ pkgs.cabextract ];
          installPhase = ''
            runHook preInstall
            cabextract -q "$src" -d outer
            cabextract -q outer/vc_red.cab -d inner
            mkdir -p "$out"
            suffix=30729.6161
            cp inner/mfc90.dll.$suffix.Microsoft_VC90_MFC_x64.QFE      "$out/mfc90.dll"
            cp inner/mfc90u.dll.$suffix.Microsoft_VC90_MFC_x64.QFE     "$out/mfc90u.dll"
            cp inner/msvcp90.dll.$suffix.Microsoft_VC90_CRT_x64.QFE    "$out/msvcp90.dll"
            cp inner/msvcr90.dll.$suffix.Microsoft_VC90_CRT_x64.QFE    "$out/msvcr90.dll"
            cp inner/manifest.$suffix.Microsoft_VC90_CRT_x64.QFE       "$out/Microsoft.VC90.CRT.manifest"
            cp inner/manifest.$suffix.Microsoft_VC90_MFC_x64.QFE       "$out/Microsoft.VC90.MFC.manifest"
            runHook postInstall
          '';
          meta.license = pkgs.lib.licenses.unfree;
        };

        # FEMM 4.2 (Windows binaries, run under Wine). Extracted natively with
        # innoextract so the package itself is reproducible without Wine.
        femm = pkgs.stdenvNoCC.mkDerivation {
          pname = "femm-bin";
          version = "4.2-21apr2019";

          src = pkgs.fetchurl {
            url = "https://www.femm.info/doku/lib/exe/fetch.php?media=upload:files:femm42bin_x64_21apr2019.exe";
            name = "femm42bin_x64_21apr2019.exe";
            sha256 = "0hr0jlpcmnni3hahrcs9drmjc6h9xyv7a7zgzz45wc5nj24lwf0p";
          };

          dontUnpack = true;
          nativeBuildInputs = [ pkgs.innoextract ];

          installPhase = ''
            runHook preInstall
            innoextract --extract --output-dir ext "$src"
            mkdir -p "$out"
            cp -r ext/app/* "$out"/
            cp ${vc90-runtime}/* "$out"/bin/
            runHook postInstall
          '';

          meta = {
            description = "Finite Element Method Magnetics (Windows binary distribution)";
            homepage = "https://www.femm.info";
            license = pkgs.lib.licenses.aladdin; # Aladdin Free Public License
            platforms = [ "x86_64-linux" ];
          };
        };

        # Headless FEMM Lua runner: femm-lua <script.lua> [more args...]
        # FEMM's app dir is copied to a writable location (.femm/) on first run
        # because FEMM writes state next to its binaries; Wine prefix lives there too.
        # FEMM is x64-only, so 64-bit-only Wine suffices (and is much lighter
        # than the wow64 build).
        femm-lua = pkgs.writeShellApplication {
          name = "femm-lua";
          runtimeInputs = [ pkgs.wine64Packages.stable pkgs.xvfb-run ];
          text = ''
            : "''${FEMM_HOME:=$PWD/.femm}"
            export WINEPREFIX="''${WINEPREFIX:-$FEMM_HOME/wineprefix}"
            export WINEDEBUG="''${WINEDEBUG:--all}"
            export WINEDLLOVERRIDES="mscoree=,mshtml=" # no mono/gecko install prompts

            if [ ! -e "$FEMM_HOME/app/bin/femm.exe" ]; then
              mkdir -p "$FEMM_HOME"
              cp -r --no-preserve=mode,ownership ${femm} "$FEMM_HOME/app"
            fi

            script="$(realpath "$1")"
            shift
            # Wine maps Z: to /; FEMM accepts forward slashes in paths.
            exec xvfb-run -a wine "$FEMM_HOME/app/bin/femm.exe" \
              "-lua-script=Z:$script" -windowhide "$@"
          '';
        };

        default = femm;
      });

      devShells = forAllSystems (pkgs: {
        default = pkgs.mkShell {
          packages = [
            self.packages.${pkgs.stdenv.hostPlatform.system}.femm-lua
            pkgs.wine64Packages.stable
            pkgs.xvfb-run
            (pkgs.python3.withPackages (ps: with ps; [
              numpy
              scipy
              pandas
              matplotlib
            ]))
          ];
          env.FEMM_DIR = "${self.packages.${pkgs.stdenv.hostPlatform.system}.femm}";
        };
      });
    };
}
