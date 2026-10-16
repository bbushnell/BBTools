package prok;

/** Stable resource identities and calibrated policies; only the 5S pilot is registered.
 * Other families retain their existing factories until individually migrated.
 * @author Brian, Raiden */
public final class NcrnaFamilyRegistry {

	public static Entry get(String id){
		if("5S".equalsIgnoreCase(id) || "prok5S".equalsIgnoreCase(id)){return PROK5S;}
		throw new IllegalArgumentException("No unified ncRNA resource registration for family: "+id);
	}

	public static final class Entry {
		private Entry(String id_, String directory_, String member_, String gate_, String[] models_){
			assert(id_!=null && directory_!=null && models_.length>0) : "Resource identity and ordered models bind every loaded family";
			id=id_;directory=directory_;member=member_;gate=gate_;models=models_;
		}
		public final String id, directory, member, gate;
		final String[] models;
		final String[] coreRoles={"consensus", "seeds"};
		final String[] endpointRoles={"tables", "endpoints"};
		final int seedK=17, minLen=90, windowPad=150, indexK=7, topN=3, minHits=2;
		final float scoreA=0f, weight=3f, idPass=.80f, idBorderline=.60f, hbmPass=.60f, collapse=.85f;
		final int outputType=ProkObject.r5S, seedMinHits=2;
	}
	private static final Entry PROK5S=new Entry("prok5S", "callgenes_5S", "", "prok5S seeded selection",
		new String[]{"universal", "bacteria", "archaea", "plastid"});
	private NcrnaFamilyRegistry(){}
}
