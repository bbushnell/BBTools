package clade;

import tax.TaxTree;

/**
 * Quick test harness for CladeConfidence NN loading and prediction.
 * Tests loading time, bin index, and confidence at various levels/lengths.
 */
public class ConfidenceTest {

	public static void main(String[] args) {

		System.out.println("=== CladeConfidence Test ===\n");

		// Time the static initialization (loads 88 NNs + calibration)
		long t0=System.nanoTime();
		float dummy=CladeConfidence.probCorrect(100000, 0.01f, 0.05f, 0.001f, 0.0005f,
				0.1f, 0.2f, 0.3f, TaxTree.SPECIES);
		long t1=System.nanoTime();
		System.out.printf("NN load + first prediction: %.1f ms%n", (t1-t0)/1e6);

		// Time subsequent predictions (should be fast)
		t0=System.nanoTime();
		int reps=10000;
		for(int i=0; i<reps; i++){
			CladeConfidence.probCorrect(100000, 0.01f, 0.05f, 0.001f, 0.0005f,
					0.1f, 0.15f, 0.2f, TaxTree.GENUS);
		}
		t1=System.nanoTime();
		System.out.printf("%d predictions: %.1f ms (%.1f us each)%n%n", reps, (t1-t0)/1e6, (t1-t0)/1e3/reps);

		int[] lengths={1000, 3000, 7000, 15000, 50000, 100000, 200000, 500000, 1000000, 2000000, 4000000};

		// Test confidence at different levels for a "good" hit (low difs)
		System.out.println("--- Good hit (low difs), len=500K ---");
		float gcdif=0.005f, strdif=0.02f, hhdif=0.001f, cagadif=0.0003f;
		float k3dif=0.04f, k4dif=0.06f, k5dif=0.08f;
		printConfidence(500000, gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif);

		// Test "medium" hit
		System.out.println("--- Medium hit, len=500K ---");
		k3dif=0.15f; k4dif=0.25f; k5dif=0.35f;
		printConfidence(500000, gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif);

		// Test "bad" hit (high difs)
		System.out.println("--- Bad hit (high difs), len=500K ---");
		k3dif=0.3f; k4dif=0.5f; k5dif=0.7f;
		gcdif=0.05f; strdif=0.15f;
		printConfidence(500000, gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif);

		// Test length effect on same difs
		System.out.println("--- Same difs (k5=0.15), varying length ---");
		gcdif=0.01f; strdif=0.03f; hhdif=0.001f; cagadif=0.0003f;
		k3dif=0.08f; k4dif=0.12f; k5dif=0.15f;
		int[] testLens={3000, 10000, 50000, 200000, 1000000, 4000000};
		for(int len : testLens){
			float sp=CladeConfidence.probCorrect(len, gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif, TaxTree.SPECIES);
			float ge=CladeConfidence.probCorrect(len, gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif, TaxTree.GENUS);
			float fa=CladeConfidence.probCorrect(len, gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif, TaxTree.FAMILY);
			float ph=CladeConfidence.probCorrect(len, gcdif, strdif, hhdif, cagadif, k3dif, k4dif, k5dif, TaxTree.PHYLUM);
			System.out.printf("  len=%,9d  sp=%.3f  ge=%.3f  fa=%.3f  ph=%.3f%n", len, sp, ge, fa, ph);
		}
		System.out.println();

		// Test sigmoid fallback (no gcdif/strdif)
		System.out.println("--- Sigmoid fallback (4-arg) ---");
		for(int len : testLens){
			float sp=CladeConfidence.probCorrect(len, k3dif, k4dif, k5dif, TaxTree.SPECIES);
			float ge=CladeConfidence.probCorrect(len, k3dif, k4dif, k5dif, TaxTree.GENUS);
			System.out.printf("  len=%,9d  sp=%.3f  ge=%.3f%n", len, sp, ge);
		}

		// Species deep-dive: what k5dif is needed for species confidence?
		System.out.println("\n--- Species confidence at various k5dif, len=3M ---");
		float[] k5vals={0.001f, 0.005f, 0.01f, 0.02f, 0.03f, 0.05f, 0.08f, 0.1f, 0.15f, 0.2f};
		for(float k5 : k5vals){
			float k4t=k5*0.8f, k3t=k5*0.5f;
			float sp=CladeConfidence.probCorrect(3000000, 0.005f, 0.02f, 0.001f, 0.0003f,
					k3t, k4t, k5, TaxTree.SPECIES);
			float ge=CladeConfidence.probCorrect(3000000, 0.005f, 0.02f, 0.001f, 0.0003f,
					k3t, k4t, k5, TaxTree.GENUS);
			System.out.printf("  k5=%.3f  sp=%.4f  ge=%.4f%n", k5, sp, ge);
		}

		// Same at 5K
		System.out.println("\n--- Species confidence at various k5dif, len=5K ---");
		for(float k5 : k5vals){
			float k4t=k5*0.8f, k3t=k5*0.5f;
			float sp=CladeConfidence.probCorrect(5000, 0.005f, 0.02f, 0.001f, 0.0003f,
					k3t, k4t, k5, TaxTree.SPECIES);
			float ge=CladeConfidence.probCorrect(5000, 0.005f, 0.02f, 0.001f, 0.0003f,
					k3t, k4t, k5, TaxTree.GENUS);
			System.out.printf("  k5=%.3f  sp=%.4f  ge=%.4f%n", k5, sp, ge);
		}

		System.out.println("\nDone.");
	}

	private static void printConfidence(int len, float gcdif, float strdif,
			float hhdif, float cagadif, float k3dif, float k4dif, float k5dif) {
		String[] names={"species","genus","family","order","class","phylum","kingdom","domain"};
		int[] levels={TaxTree.SPECIES, TaxTree.GENUS, TaxTree.FAMILY, TaxTree.ORDER,
				TaxTree.CLASS, TaxTree.PHYLUM, TaxTree.KINGDOM, TaxTree.DOMAIN};
		for(int i=0; i<levels.length; i++){
			float c=CladeConfidence.probCorrect(len, gcdif, strdif, hhdif, cagadif,
					k3dif, k4dif, k5dif, levels[i]);
			System.out.printf("  %-8s: %.4f (%.1f%%)%n", names[i], c, c*100);
		}
		System.out.println();
	}
}
