package prok;

/** Explicit experimental mitochondrial additions to the existing16S library.
 * The original five models retain their measured floors and endpoint policy.
 * Added models use raw ends and a full-model-relative support floor; that floor
 * is an experimental policy, not a newly calibrated biological cutoff.
 * @author Brian Bushnell, Raiden */
public final class MitoSsuModelPolicy {
	/** Coverage-selected candidates retain the exact original prefix; the suffix
	 * is explicitly named, not an arbitrary extension of a production library. */
	public static boolean coverage(final String family, final String[] names){
		if(!"16S".equals(family) || names==null){return false;}
		final String[] original=NcrnaFamilyRegistry.get("16S").models;
		if(names.length<=original.length){return false;}
		for(int i=0; i<names.length; i++){
			if(names[i]==null){return false;}
			final String token=RrnaResourceIO.first(names[i]);
			if(i<original.length ? !token.equals(original[i]) : !token.matches("mitochondria_coverage_[A-Za-z0-9_]+")){return false;}
		}
		return true;
	}
	public static boolean extended(final String family, final String[] names){return expanded(family,names) || coverage(family,names);}
	public static boolean expanded(final String family, final String[] names){
		if(!"16S".equals(family) || names==null){return false;}
		final String[] original=NcrnaFamilyRegistry.get("16S").models;
		if(names.length!=original.length+ADDITIONS.length){return false;}
		for(int i=0; i<names.length; i++){
			if(names[i]==null || !name(i, original).equals(RrnaResourceIO.first(names[i]))){return false;}
		}
		return true;
	}
	static String name(final int index, final String[] original){
		assert(index>=0 && index<original.length+ADDITIONS.length):"Extended model order must retain the original roster followed by the five named mitochondrial consensuses";
		return index<original.length ? original[index] : ADDITIONS[index-original.length];
	}
	static String[] rawNames(final NcrnaFamily family, final String[] original){
		if(!extended(family.name, family.modelNames)){return original;}
		final int first=NcrnaFamilyRegistry.get("16S").models.length;
		final String[] result=java.util.Arrays.copyOf(original, original.length+family.modelNames.length-first);
		for(int i=first; i<family.modelNames.length; i++){result[original.length+i-first]=RrnaResourceIO.first(family.modelNames[i]);}
		return result;
	}
	/** Full Quantum identity uses M/(M+S+I+D), with M+S+I=qLen
	 * (idaligner.Tracer.postprocess). At least pass*qLen reference bases are
	 * needed. Use the same float comparison at the integer boundary, avoiding
	 * ceil(float*length) rejecting e.g.33/50 at .66. Edge-clipped candidates
	 * deliberately retain this full-model-relative support floor too. */
	public static int minimumLength(final int modelLength, final float pass){
		if(modelLength<17 || !Float.isFinite(pass) || pass<=0f || pass>1f){
			throw new IllegalArgumentException("Mitochondrial support requires a complete model of at least17bases and identity in(0,1]: "+modelLength+", "+pass);
		}
		int n=Math.max(1, (int)(pass*modelLength));
		while(n>1 && (float)(n-1)/modelLength>=pass){n--;}
		while((float)n/modelLength<pass){n++;}
		assert(n<=modelLength && (n==1 || (float)(n-1)/modelLength<pass)):
			"Minimum support must be the first integral span meeting the full-model-relative identity fraction";
		return n;
	}
	public static final String[] ADDITIONS={"mitochondria_chordate", "mitochondria_arthropod", "mitochondria_land_plant",
		"mitochondria_insect", "mitochondria_fallback"};
	private MitoSsuModelPolicy(){}
}
