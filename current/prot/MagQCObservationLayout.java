package prot;

/**
 * Fixed observation types and item order for MAG-QC subnet inputs.
 * The five historical ncRNA items retain their offsets; structural anticodon
 * counts and the unknown residual follow them in the extended expected-copy table.
 * @author Yoimiya
 */
public final class MagQCObservationLayout {

	private MagQCObservationLayout(){}

	/** Returns the fixed special width, or -1 for a family-rank subset. */
	public static int specialCount(String type){
		if("famset".equals(type)){return -1;}
		if("ncrna".equals(type)){return LEGACY_COUNT;}
		if("rrna".equals(type)){return RRNA_COUNT;}
		if("trna_anticodon".equals(type)){return ANTICODON_COUNT;}
		throw new IllegalArgumentException("Unknown MAG-QC observation type: "+type);
	}

	/** Canonical N-item key at a non-protein offset, independent of family count. */
	public static String nonproteinKey(int index){
		if(index<0 || index>=EXTENDED_COUNT){
			throw new IllegalArgumentException("Non-protein observation offset outside [0,70): "+index);
		}
		if(index<LEGACY_COUNT){return LEGACY_KEYS[index];}
		if(index==EXTENDED_COUNT-1){return "trna_anticodon_unknown";}
		assert(index-LEGACY_COUNT<64) : "The 64 structural anticodon codes precede the unknown residual";
		return "trna_anticodon_"+(index-LEGACY_COUNT);
	}

	/** Serialized description checked independently of network width and payload hash. */
	public static String definition(String type){
		final int count=specialCount(type);
		if(count<0){return "-";}
		if("ncrna".equals(type)){return "5 ordered ncRNA observations";}
		if("rrna".equals(type)){return "r16,r23,r5,rother";}
		assert("trna_anticodon".equals(type)) : "specialCount accepts only the four frozen subnet types";
		return "structural anticodon codes 0..63,unknown residual";
	}

	public static final int LEGACY_COUNT=5;
	public static final int RRNA_COUNT=4;
	public static final int ANTICODON_COUNT=65;
	public static final int EXTENDED_COUNT=LEGACY_COUNT+ANTICODON_COUNT;
	private static final String[] LEGACY_KEYS={"r16","r23","r5","rother","trna"};
}
