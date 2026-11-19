package ddl;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import dna.Data;
import shared.Resources;

/** Exact-geometry official blacklists, shared by producers and consumers.
 * Brian selected reduced12 K6 shadow4096, reduced16 K5 raw-shadow1024 and d120 K6 nonoverlap shadow2048.
 * Other geometries retain their explicitly configured filtering behavior.
 * @author Yaoyao
 */
public final class ProteinBlacklistPolicy{

	/** Static policy lookup and loading only. */
	private ProteinBlacklistPolicy(){}

	/*--------------------------------------------------------------*/
	/*----------------       Policy Selection       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns the exact geometry's descriptor, or null for null/unrecognized geometry.
	 * Lookup does not resolve resource paths or load a filter.
	 */
	public static Descriptor policy(ProteinKmerGeometry geometry){
		if(geometry==null){return null;}
		if(REDUCED12.geometry.equals(geometry)){return REDUCED12;}
		if(REDUCED16.geometry.equals(geometry)){return REDUCED16;}
		if(DIMER120.geometry.equals(geometry)){return DIMER120;}
		return null;
	}

	/** Only the exact accepted mapping and word length select a policy. */
	public static boolean applies(ProteinKmerGeometry geometry){return policy(geometry)!=null;}

	/**
	 * Resolves an omitted official filter or validates an explicit filter's policy.
	 * Unknown geometries retain supplied after geometry validation. For DIMER120,
	 * a supplied filter is checked against descriptor metadata and returned directly.
	 * The monomer policies instead compare supplied keys with a loaded official filter;
	 * filters whose identity was not measured select the structural official loader there.
	 * @return Supplied, an official singleton, or null when no policy/filter applies
	 */
	public static DDLPentamerBlacklist resolve(ProteinKmerGeometry geometry, DDLPentamerBlacklist supplied){
		final Descriptor selected=policy(geometry);
		if(selected==null){
			if(supplied!=null){supplied.validate(geometry);}
			return supplied;
		}
		if(selected==DIMER120 && supplied!=null){
			selected.requireSupplied(supplied);
			return supplied;
		}
		if(selected.isCachedOfficial(supplied)){return supplied;}
		if(supplied!=null && !supplied.identityMeasured()){
			selected.requireSupplied(supplied);
			final DDLPentamerBlacklist official=selected.structuralOfficial();
			supplied.requireOfficial(official, geometry, selected.label);
			return official;
		}
		final DDLPentamerBlacklist official=selected.official();
		if(supplied!=null && supplied!=official){supplied.requireOfficial(official, geometry, selected.label);}
		return official;
	}

	/** Explicit experimental policies never select a default by class count.
	 * Structural identity is provenance only; experiment controllers bind exact files.
	 * Requires an absolute path, structural=true and complete policy metadata; an
	 * exact official geometry must retain that descriptor's identity/count/mode/provenance.
	 */
	public static DDLPentamerBlacklist loadExplicit(ProteinKmerGeometry geometry, String path, String identity,
			int words, int filterK, String modeText, String provenance, boolean structural)throws Exception{
		if(geometry==null || path==null || identity==null || words<1 || filterK<1 || modeText==null ||
			provenance==null || provenance.trim().isEmpty() || !structural){
			throw new IllegalArgumentException("Explicit clustering blacklist requires structural verification, path, identity, words, k, mode and provenance");
		}
		final DDLPentamerBlacklist.MatchMode mode=DDLPentamerBlacklist.MatchMode.parse(modeText);
		final Descriptor official=policy(geometry);
		if(official!=null && (words!=official.words || filterK!=official.filterK || mode!=official.mode ||
			!identity.equals(official.sha80) || !provenance.equals(official.provenance))){
			throw new IllegalArgumentException("Explicit policy cannot replace the selected "+official.label);
		}
		if(mode==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL && filterK!=geometry.k()){
			throw new IllegalArgumentException("Reduced blacklist word length must match the complete sketch geometry");
		}
		final Path input=Paths.get(path);
		if(!input.isAbsolute()){throw new IllegalArgumentException("Explicit structural blacklist requires an absolute path");}
		final DDLPentamerBlacklist filter=mode==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL ?
			DDLPentamerBlacklist.loadReducedStructural(input, identity, words, geometry, provenance) :
			DDLPentamerBlacklist.loadStructural(input, identity, words, filterK, mode, provenance);
		filter.validate(geometry);
		return resolve(geometry, filter);
	}

	/** Requires the official identity token when geometry selects a policy.
	 * This checks metadata only; it does not load a file or retrofit stored sketches.
	 */
	public static void requireStored(ProteinKmerGeometry geometry, String identity){
		final Descriptor selected=policy(geometry);
		if(selected!=null && !selected.sha80.equals(identity)){
			throw new IllegalArgumentException(selected.label+" blacklist is required; rebuild/recluster historical models with the same filter");
		}
	}

	/** Compatibility entry point: this still resolves the original 12K6 resource. */
	public static Path resource(){return REDUCED12.resource();}

	/*--------------------------------------------------------------*/
	/*----------------       Policy Descriptor      ----------------*/
	/*--------------------------------------------------------------*/

	/** One immutable description supplies the geometry, matching space and persisted metadata. */
	public static final class Descriptor{
		/** Binds one exact mapping and its expected resource metadata; performs no I/O. */
		private Descriptor(ProteinKmerGeometry geometry_, int words_, DDLPentamerBlacklist.MatchMode mode_,
				String sha80_, String resourceName_, String provenance_, String label_){
			this(geometry_, words_, mode_, sha80_, sha80_, resourceName_, resourceName_, provenance_, label_);
		}
		/** Binds one exact mapping, separating shipped-file bytes from persisted policy identity. */
		private Descriptor(ProteinKmerGeometry geometry_, int words_, DDLPentamerBlacklist.MatchMode mode_,
				String sha80_, String fileSha80_, String resourceName_, String provenance_, String label_){
			this(geometry_, words_, mode_, sha80_, fileSha80_, resourceName_, "bbgap/"+resourceName_, provenance_, label_);
		}
		/** Binds one exact mapping and its installed resource lookup path. */
		private Descriptor(ProteinKmerGeometry geometry_, int words_, DDLPentamerBlacklist.MatchMode mode_,
				String sha80_, String fileSha80_, String resourceName_, String resourcePath_, String provenance_, String label_){
			geometry=geometry_;
			filterK=geometry.k();
			words=words_;
			mode=mode_;
			sha80=sha80_;
			fileSha80=fileSha80_;
			resourceName=resourceName_;
			resourcePath=resourcePath_;
			provenance=provenance_;
			label=label_;
		}

		/** Resolve a portable checkout or installed BBTools resource. */
		public Path resource(){
			final String path=Data.findPath("?"+resourcePath, false);
			if(path==null){
				final ArrayList<String> missing=new ArrayList<String>(1);
				missing.add("?"+resourcePath);
				Resources.throwMissing(missing);
			}
			return Paths.get(path).toAbsolutePath().normalize();
		}

		/** Loads this policy's format, projection and exact word count.
		 * Structural mode retains the declared identity without measuring a digest;
		 * false selects the legacy exact-pin loader. The monomer resolve() branch
		 * additionally compares supplied keys with its official singleton; the
		 * supplied-DIMER120 branch checks descriptor metadata instead.
		 */
		public DDLPentamerBlacklist load(Path path, boolean structural)throws Exception{
			final DDLPentamerBlacklist filter;
			if(mode==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL){
				filter=structural ? DDLPentamerBlacklist.loadReducedStructural(path, sha80, words, geometry, provenance) :
					loadPublished(path);
			}else{
				filter=structural ? DDLPentamerBlacklist.loadStructural(path, sha80, words, filterK, mode, provenance) :
					loadPublished(path);
			}
			if(filter.size()!=words){throw new IllegalArgumentException(label+" blacklist cardinality differs from the selected policy");}
			filter.validate(geometry);
			return filter;
		}

		/** Returns the lazily initialized legacy pinned singleton for this exact descriptor. */
		private DDLPentamerBlacklist official(){
			DDLPentamerBlacklist filter=officialFilter;
			if(filter!=null){return filter;}
			synchronized(this){
				filter=officialFilter;
				if(filter==null){officialFilter=filter=loadOfficial();}
			}
			return filter;
		}
		/** Structural callers compare all canonical keys without promoting inherited pins to measured identity. */
		private DDLPentamerBlacklist structuralOfficial(){
			DDLPentamerBlacklist filter=structuralFilter;
			if(filter!=null){return filter;}
			synchronized(this){
				filter=structuralFilter;
				if(filter==null){structuralFilter=filter=loadStructuralOfficial();}
			}
			return filter;
		}
		private boolean isCachedOfficial(DDLPentamerBlacklist filter){return filter!=null && filter==officialFilter;}
		/** Checks the supplied filter's descriptor metadata and geometry, without loading a reference. */
		private void requireSupplied(DDLPentamerBlacklist supplied){
			if(!sha80.equals(supplied.sha80()) || supplied.size()!=words || supplied.wordLength()!=filterK ||
					supplied.matchMode()!=mode || !provenance.equals(supplied.provenance())){
				throw new IllegalArgumentException(label+" supplied blacklist differs from the selected policy");
			}
			supplied.validate(geometry);
		}
		/** Loads the legacy pinned singleton, wrapping loader failures with policy context. */
		private DDLPentamerBlacklist loadOfficial(){
			try{return load(resource(), false);}
			catch(Resources.MissingResourceException e){throw e;}
			catch(Exception e){throw new IllegalArgumentException("Could not load the mandatory "+label+" blacklist", e);}
		}
		/** Loads the installed resource structurally for comparing inherited identities. */
		private DDLPentamerBlacklist loadStructuralOfficial(){
			try{return load(resource(), true);}
			catch(Resources.MissingResourceException e){throw e;}
			catch(Exception e){throw new IllegalArgumentException("Could not load the mandatory structural "+label+" blacklist", e);}
		}
		/** Verifies the published resource bytes while retaining the historical policy identity. */
		private DDLPentamerBlacklist loadPublished(Path path)throws Exception{
			return DDLPentamerBlacklist.loadPublished(path, fileSha80, sha80, words, filterK, mode, provenance,
				mode==DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL ? geometry : null);
		}
		/** Exact word-projection geometry selected by this policy. */
		public final ProteinKmerGeometry geometry;
		/** Required word length and exact number of distinct filter words. */
		public final int filterK, words;
		/** Whether matching uses raw literals, prefixes or projected full words. */
		public final DDLPentamerBlacklist.MatchMode mode;
		/** Exact shipped file bytes; can differ from the pre-header policy identity. */
		public final String fileSha80;
		/** Declared identity, portable resource filename, provenance token and diagnostic label. */
		public final String sha80, resourceName, provenance, label;
		/** Installed BBTools resource path. */
		private final String resourcePath;
		private volatile DDLPentamerBlacklist officialFilter, structuralFilter;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Fixed Policies        ----------------*/
	/*--------------------------------------------------------------*/

	// Existing public constants retain their original 12K6 meaning.
	public static final String ALPHABET="all:A/C/DE/FWY/G/IV/LM/NH/P/Q/RK/ST";
	public static final ProteinKmerGeometry GEOMETRY=ProteinKmerGeometry.make(6,ALPHABET);
	public static final int WORDS=4096;
	public static final String SHA80="3a2ac44cbac771f86fd1";
	public static final String RESOURCE="protein_blacklist_reduced12_k6_shadow4096.faa";
	public static final String PROVENANCE="official_reduced12_k6_shadow4096;selector=20b7e537d848652ec427c872ad341ef5bf1212ad";
	private static final Descriptor REDUCED12=new Descriptor(GEOMETRY,WORDS,DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL,
		SHA80,RESOURCE,PROVENANCE,"12K6 official reduced-shadow4096");
	private static final Descriptor REDUCED16=new Descriptor(
		ProteinKmerGeometry.make(5,"all:A/C/D/E/FWY/G/H/IV/K/LM/N/P/Q/R/S/T"),1024,DDLPentamerBlacklist.MatchMode.LITERAL,
		"f731f9dfe82458b93234","0e368b6d191d317d7bfe","protein_blacklist_reduced16_k5_raw_shadow1024.faa.gz",
		"official_reduced16_k5_raw_shadow1024;selector=c77644d0fc8960d2289292a248bd1f32c8270ed9","16K5 official raw-shadow1024");
	/** Exact ordered mapping from the accepted d120 K6 nonoverlap candidate recipe. */
	public static final String DIMER120_ALPHABET="dimer:1,3,5:AA/AC,RC,NC,DC,EC,QC,GC,HC,KC,PC,SC,TC/AD,SD,TD/AE,AQ/AF,RF,KF,SF,TF/AG,SG,TG/"+
		"AH,RH,EH,QH,KH,SH,TH/AI,AV/AK,SK,TK/AL/AM,RM,QM,KM,SM,TM/AN/AP,RP,EP,QP,KP,SP,TP/AR/AS,AT/"+
		"AW,RW,NW,DW,EW,QW,HW,KW,SW,TW/AY,RY,DY,EY,KY,SY,TY/CA,CR,CN,CD,CE,CQ,CG,CH,CK,CP,CS,CT/"+
		"CC,IC,LC,MC,FC,WC,YC,VC/CI,CL,CM,CF,CW,CY,CV/DA,DS,DT/DD,ED/DE,EE/DF,EF/DI,DL,DV/DM,EM/DQ,EQ/DR,DK/"+
		"EA,QA/EI,EL,EV/ER,EK/ES,ET/FA,FR,FE,FQ,FK,FS,FT/FF,FY,YF,YY/FG,YG/FI,FL,FV/FN,FD,YN,YD/FP,WP,YP/"+
		"GA,GS,GT/GF,GW,GY/GG/GI,GL,GM,GV/GN,GD/GP/GR,GE,GQ,GK/HA,HR,HE,HQ,HK,HS,HT/HI,HL,HM,HF,HY,HV/"+
		"HN,HD,HH/IA,VA/ID,LD,VD/IE,LE,VE/IF,LF,VF/IG,LG,MG,VG/IH,LH,MH,FH,WH,YH,VH/II,IV,VI,VV/IK,VK/"+
		"IL,LL,VL/IM,LM,MM,FM,YM,VM/IN,LN,VN/IP,LP,MP,VP/IQ,LQ,VQ/IR,VR/IS,VS/IT,VT/IW,LW,MW,FW,WW,YW,VW/"+
		"IY,LY,VY/KA,KS,KT/KD,KE/KI,KL,KV/LA/LI,LV/LR,LK/LS,LT/MA,MR,ME,MQ,MK,MS,MT/MI,ML,MF,MY,MV/MN,MD/"+
		"NA,NS,NT/ND,NE/NF,NY/NG,DG,EG,QG/NH,DH,GH/NI,NL,NM,NV/NN,DN,EN,SN,TN/NP,DP,HP/NR,NQ,NK/"+
		"PA,PS,PT/PD,PE/PF,PW,PY/PG/PI,PL,PM,PV/PP/PR,PN,PQ,PH,PK/QD,QE/QF,QY/QI,QL,QV/QR,QQ,QK/QS,QT/"+
		"RA/RD,RE/RG,HG,KG/RI,RL,RV/RK,KK/RN,QN,KN/RQ,KQ/RR,KR/RS,RT/SA/SE,TE/SI,SV/SL,TL/SQ,TQ/SR,TR/"+
		"SS,TS/ST,TT/TA/TI,TV/WA,WR,WN,WD,WE,WQ,WG,WK,WS,WT/WI,WL,WM,WF,WY,WV/YA,YR,YE,YQ,YK,YS,YT/YI,YL,YV";
	private static final Descriptor DIMER120=new Descriptor(ProteinKmerGeometry.make(6,DIMER120_ALPHABET),2048,
		DDLPentamerBlacklist.MatchMode.REDUCED_LITERAL,"c0fb8ec0a91bb50535c5","73c76063bac23e79aa33",
		"protein_blacklist_d120_k6_nonoverlap_shadow2048.faa.gz",
		"official_d120_k6_nonoverlap_shadow2048;selector=8126cccd270e71d3ead11e2d7b2f6416702bad0f",
		"d120 K6 nonoverlap official reduced-shadow2048");
}
