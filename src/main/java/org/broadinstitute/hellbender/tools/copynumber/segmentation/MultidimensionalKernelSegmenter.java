package org.broadinstitute.hellbender.tools.copynumber.segmentation;

import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.List;

/**
 * Segments copy-ratio and alternate-allele-fraction data using kernel segmentation.  Segments do not span chromosomes.
 * Only the first allele-fraction site in each copy-ratio interval is used.  The alternate-allele fraction in
 * copy-ratio intervals that do not contain any sites is imputed to be balanced at 0.5.
 * Refactored to be a thin wrapper around the {@link MultisampleMultidimensionalKernelSegmenter}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class MultidimensionalKernelSegmenter {
    private final MultisampleMultidimensionalKernelSegmenter segmenter;

    public MultidimensionalKernelSegmenter(final CopyRatioCollection denoisedCopyRatios,
                                           final AllelicCountCollection allelicCounts) {
        Utils.nonNull(denoisedCopyRatios);
        Utils.nonNull(allelicCounts);
        segmenter = new MultisampleMultidimensionalKernelSegmenter(
                Collections.singletonList(denoisedCopyRatios),
                Collections.singletonList(allelicCounts));
    }

    /**
     * Segments the internally held {@link CopyRatioCollection} and {@link AllelicCountCollection}
     * using a separate {@link KernelSegmenter} for each chromosome.
     * @param kernelVarianceCopyRatio       variance of the Gaussian kernel used for copy-ratio data;
     *                                      if zero, a linear kernel is used instead
     * @param kernelVarianceAlleleFraction  variance of the Gaussian kernel used for allele-fraction data;
     *                                      if zero, a linear kernel is used instead
     * @param kernelScalingAlleleFraction   relative scaling S of the kernel K_AF for allele-fraction data
     *                                      to the kernel K_CR for copy-ratio data;
     *                                      the total kernel is K_CR + S * K_AF
     */
    public SimpleIntervalCollection findSegmentation(final int maxNumSegmentsPerChromosome,
                                                     final double kernelVarianceCopyRatio,
                                                     final double kernelVarianceAlleleFraction,
                                                     final double kernelScalingAlleleFraction,
                                                     final int kernelApproximationDimension,
                                                     final List<Integer> windowSizes,
                                                     final double numChangepointsPenaltyLinearFactor,
                                                     final double numChangepointsPenaltyLogLinearFactor) {
        return segmenter.findSegmentation(
                maxNumSegmentsPerChromosome,
                kernelVarianceCopyRatio,
                kernelVarianceAlleleFraction,
                kernelScalingAlleleFraction,
                kernelApproximationDimension,
                windowSizes,
                numChangepointsPenaltyLinearFactor,
                numChangepointsPenaltyLogLinearFactor);
    }
}
