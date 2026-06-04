# M33-Corbelli Remaining Limitations Before Physics Execution

**Lane label:** M33-Corbelli source-native thick-disk reconstruction  
**Status:** Format B reconstructed board; mass closure passed; overlay benchmark passed; annular thick-disk kernel provisionally frozen; final board-freeze package complete; age status Amber with bracketed $t_{50}$ locked for sensitivity use; no Hermes/MOND execution authorized.

## Remaining limitations

1. **Source-native, not SPARC-equivalent.** The stellar mass profile is Corbelli's SPS / pixel-SED $\Sigma_\star(R)$ with radially varying M/L provenance, not SPARC unit-M/L 3.6 μm photometry.

2. **Component velocities are reconstructed products.** $V_{\rm gas}(R)$ and $V_{\rm disk}(R)$ are computed from public surface-density profiles using the provisionally frozen annular thick-disk kernel. They are not direct source-published rotmod arrays.

3. **Kernel is provisionally frozen for this lane only.** The annular quadrature thick-disk kernel passed the source Figure 12 top-panel overlay benchmark and is accepted for the M33-Corbelli source-native lane. It should not be generalized automatically to other lanes.

4. **H I profile mass excess is accepted but documented.** The integrated H I profile mass exceeds the true H I mass by about 12%, consistent with the source warning about outer-ring filling. This is not treated as a failure, but it remains part of the provenance record.

5. **The board is 1D/axisymmetric.** The reconstruction does not implement a full 3D warp model or non-axisymmetric structure. It follows the source-native radial mass-model treatment.

6. **Age status is Amber, bracketed.** No direct global published M33 $t_{50}$ was identified. The locked treatment is $t_{50}=6.0\pm1.2$ Gyr primary/source-native and $t_{50}=7.1\pm1.3$ Gyr conservative older-weighted. Both must be reported in any eventual Hermes result.

7. **No model result exists.** No Hermes run, MOND run, or gravity comparison has been performed or authorized.
