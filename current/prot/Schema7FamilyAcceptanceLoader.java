package prot;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import fileIO.ByteFile;
import parse.LineParser1;

/**
 * Strict runtime loader for the sealed schema-7 family-threshold artifact.
 *
 * <p>The loader first runs the maintained independent schema verifier, then
 * binds every live production resource named by the artifact.  It parses only
 * the effective threshold columns used by assignment; empirical columns and
 * floor algebra remain guarded by {@link Schema7ThresholdVerifier}.</p>
 *
 * @author Yoimiya
 */
final class Schema7FamilyAcceptanceLoader {

	private Schema7FamilyAcceptanceLoader(){}

	/** Loads, verifies, and returns an immutable-by-ownership schema-7 table. */
	static Loaded load(final String artifactPath, final String expectedArtifactSha80,
			final String rosterPath,
			final String consensusPath, final String roleManifestPath,
			final String coreCoordinatesPath, final String coveringSetsPath,
			final String sidecarPath, final String hbmBundlePath,
			final String hbmProvenancePath){
		final String artifactSha80=DigestSuffix.requireSuffix(expectedArtifactSha80,
			"schema-7 artifact sha80");
		checkHash(artifactPath,artifactSha80,"threshold artifact");
		final Parsed parsed=parse(artifactPath);
		verifyArtifact(artifactPath,artifactSha80,parsed.header);
		final int families=positiveInt(parsed.header,"families",artifactPath);
		final long members=positiveLong(parsed.header,"members",artifactPath);
		if(parsed.rows.size()!=families){
			throw new IllegalArgumentException("Schema-7 row count "+parsed.rows.size()+
				" != header families "+families+": "+artifactPath);
		}

		final String rosterSha80=hash(parsed.header,"roster_sha80",artifactPath);
		checkHash(rosterPath,rosterSha80,"roster");
		final RuntimeRoster roster=readRoster(rosterPath,rosterSha80,families);

		final String[] repIds=new String[families];
		final String[] calibrationStatus=new String[families];
		final int[] familyIds=new int[families],n=new int[families],lengthLo=new int[families],
			lengthHi=new int[families],minRaw=new int[families],minKmer=new int[families];
		final boolean[] enabled=new boolean[families];
		final double[] minIdentity=new double[families],minR=new double[families],
			minMutual=new double[families],minKmerDensity=new double[families];
		final float[] minHbmPath=new float[families];
		long observedMembers=0;
		for(int rank=0; rank<families; rank++){
			final String[] row=parsed.rows.get(rank);
			final int observedRank=canonicalInt(row[0],"rank",rank,artifactPath);
			familyIds[rank]=canonicalInt(row[1],"family_id",rank,artifactPath);
			repIds[rank]=row[2];
			n[rank]=canonicalInt(row[3],"n",rank,artifactPath);
			calibrationStatus[rank]=row[4];
			enabled[rank]=canonicalBoolean(row[5],"acceptance_enabled",rank,artifactPath);
			lengthLo[rank]=canonicalInt(row[19],"length_lo_effective",rank,artifactPath);
			lengthHi[rank]=canonicalInt(row[22],"length_hi_effective",rank,artifactPath);
			minRaw[rank]=canonicalInt(row[7],"raw_score_effective",rank,artifactPath);
			minIdentity[rank]=canonicalFloat(row[10],"identity_effective",rank,artifactPath);
			minR[rank]=canonicalFloat(row[13],"R_effective",rank,artifactPath);
			minMutual[rank]=canonicalFloat(row[16],"mutual_overlap_effective",rank,artifactPath);
			minKmer[rank]=canonicalInt(row[25],"kmer_count_effective",rank,artifactPath);
			minKmerDensity[rank]=canonicalDouble(row[28],"kmer_density_effective",rank,artifactPath);
			minHbmPath[rank]=canonicalFloat(row[31],"hbm_path_effective",rank,artifactPath);
			if(observedRank!=rank || familyIds[rank]!=roster.familyIds[rank] ||
					!repIds[rank].equals(roster.repIds[rank]) ||
					n[rank]!=roster.members[rank] ||
					!calibrationStatus[rank].equals(roster.status[rank])){
				throw new IllegalArgumentException("Schema-7 threshold/roster mismatch at rank "+rank+".");
			}
			observedMembers=Math.addExact(observedMembers,n[rank]);
		}
		if(observedMembers!=members){
			throw new IllegalArgumentException("Schema-7 row member sum "+observedMembers+
				" != header members "+members+".");
		}

		final String roleSha80=hash(parsed.header,"role_manifest_sha80",artifactPath);
		checkHash(roleManifestPath,roleSha80,"role manifest");
		final FamilyRoleManifest roles=FamilyRoleManifest.loadSchema7(roleManifestPath,
			familyIds,repIds,rosterSha80,members);
		checkHash(roleManifestPath,roleSha80,"role manifest changed");

		final String consensusSha80=hash(parsed.header,"consensus_sha80",artifactPath);
		final String coreSha80=hash(parsed.header,"core_sha80",artifactPath);
		final String coveringSetsSha80=hash(parsed.header,"covering_sets_sha80",artifactPath);
		final String sidecarSha80=hash(parsed.header,"sidecar_sha80",artifactPath);
		final String hbmBundleSha80=hash(parsed.header,"hbm_bundle_sha80",artifactPath);
		final String hbmProvenanceSha80=hash(parsed.header,"hbm_provenance_sha80",artifactPath);
		checkHash(consensusPath,consensusSha80,"consensus");
		checkHash(coreCoordinatesPath,coreSha80,"core coordinates");
		checkHash(coveringSetsPath,coveringSetsSha80,"covering sets");
		checkHash(sidecarPath,sidecarSha80,"shortlist sidecar");
		checkStoredHash(hbmBundlePath,hbmBundleSha80,"HBM bundle");
		checkHash(hbmProvenancePath,hbmProvenanceSha80,"HBM provenance");

		final FamilyAcceptanceConfig.SidecarHeaderFields sidecar=
			FamilyAcceptanceConfig.readSidecarHeaderFields(sidecarPath);
		if(!consensusSha80.equals(sidecar.consensusSha256) ||
				!coveringSetsSha80.equals(sidecar.kmerSetsSha256)){
			throw new IllegalArgumentException("Schema-7 sidecar resource bindings do not match the sealed consensus/covering sets.");
		}
		final byte[][] hbmProvenance=HbmBundleLoader.loadSemanticProvenance(hbmProvenancePath);
		final String hbmFamilyListSha80=hbmFamilyListSha80(repIds,n);
		if(!hbmFamilyListSha80.equals(DigestSuffix.fromDigest(hbmProvenance[10])) ||
				!consensusSha80.equals(DigestSuffix.fromDigest(hbmProvenance[11]))){
			throw new IllegalArgumentException("Schema-7 HBM provenance does not bind the "+
				"sealed family-list projection and consensus.");
		}

		checkHash(rosterPath,rosterSha80,"roster changed");
		checkHash(consensusPath,consensusSha80,"consensus changed");
		checkHash(coreCoordinatesPath,coreSha80,"core coordinates changed");
		checkHash(coveringSetsPath,coveringSetsSha80,"covering sets changed");
		checkHash(sidecarPath,sidecarSha80,"shortlist sidecar changed");
		checkStoredHash(hbmBundlePath,hbmBundleSha80,"HBM bundle changed");
		checkHash(hbmProvenancePath,hbmProvenanceSha80,"HBM provenance changed");
		checkHash(artifactPath,artifactSha80,"threshold artifact changed");

		return new Loaded(artifactSha80,hash(parsed.header,"threshold_table_sha80",artifactPath),
			rosterSha80,roster.coreFamilyListSha80,consensusSha80,roleSha80,
			coreSha80,coveringSetsSha80,
			sidecarSha80,hbmBundleSha80,hbmProvenanceSha80,sidecar.coveringAlphabet,
			sidecar.coveringK,members,familyIds,repIds,n,calibrationStatus,enabled,
			lengthLo,lengthHi,minRaw,minIdentity,minR,minMutual,minKmer,
			minKmerDensity,minHbmPath,roles);
	}

	/** Reads the exact roster-v4 columns used to bind schema-7 rows. */
	private static RuntimeRoster readRoster(final String path, final String rosterSha80,
			final int families){
		final int[] familyIds=new int[families],members=new int[families];
		final String[] repIds=new String[families],status=new String[families];
		final ByteFile input=ByteFile.makeByteFile(path,false);
		final LineParser1 parser=new LineParser1((byte)'\t');
		String rosterSchema=null,supersededRosterSha80=null;
		boolean columns=false;
		int rank=0;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0){throw new IllegalArgumentException("Blank schema-7 roster line: "+path);}
				if(line[0]=='#'){
					if(columns){throw new IllegalArgumentException("Late schema-7 roster metadata: "+path);}
					final String metadata=new String(line,StandardCharsets.US_ASCII);
					final int tab=metadata.indexOf('\t');
					if(tab>1){
						final String key=metadata.substring(0,tab),value=metadata.substring(tab+1);
						if("#schema".equals(key)){
							if(rosterSchema!=null){throw new IllegalArgumentException("Duplicate roster schema: "+path);}
							rosterSchema=value;
						}else if("#supersedes_roster_v3_sha80".equals(key)){
							if(supersededRosterSha80!=null){
								throw new IllegalArgumentException("Duplicate superseded-roster identity: "+path);
							}
							supersededRosterSha80=DigestSuffix.requireSuffix(value,
								"supersedes_roster_v3_sha80");
						}
					}
					continue;
				}
				if(!columns){
					final String observed=new String(line,StandardCharsets.US_ASCII);
					if(!ROSTER_COLUMNS.equals(observed)){throw new IllegalArgumentException("Unexpected schema-7 roster columns: "+observed);}
					columns=true; continue;
				}
				if(rank>=families){throw new IllegalArgumentException("Schema-7 roster has too many rows: "+path);}
				parser.set(line);
				if(parser.terms()!=18){throw new IllegalArgumentException("Schema-7 roster row width "+parser.terms()+" != 18 at rank "+rank);}
				final int observedRank=canonicalInt(parser.parseString(0),"active_index",rank,path);
				familyIds[rank]=canonicalInt(parser.parseString(1),"family_id",rank,path);
				repIds[rank]=parser.parseString(2);
				members[rank]=canonicalInt(parser.parseString(6),"member_count",rank,path);
				final int retained=canonicalInt(parser.parseString(16),"retained_n",rank,path);
				status[rank]=parser.parseString(17);
				if(observedRank!=rank || familyIds[rank]<0 || repIds[rank].length()==0 ||
						members[rank]<1 || members[rank]!=retained){
					throw new IllegalArgumentException("Invalid schema-7 roster row at rank "+rank+".");
				}
				rank++;
			}
		}finally{input.close();}
		if(!columns || rank!=families){throw new IllegalArgumentException("Schema-7 roster row count "+rank+" != "+families+": "+path);}
		if(!ROSTER_SCHEMA.equals(rosterSchema)){
			throw new IllegalArgumentException("Schema-7 requires roster schema '"+
				ROSTER_SCHEMA+"', observed '"+rosterSchema+"': "+path);
		}
		if(supersededRosterSha80==null){
			throw new IllegalArgumentException("Schema-7 roster is missing #supersedes_roster_v3_sha80: "+path);
		}
		return new RuntimeRoster(familyIds,repIds,members,status,supersededRosterSha80);
	}

	/**
	 * Reconstructs the canonical HBM family-list input from the schema-7 roster
	 * fields that it represents. The HBM provenance binds this three-column
	 * projection, while schema 7 separately pins the complete roster-v4 file.
	 */
	static String hbmFamilyListSha80(final String[] repIds, final int[] occurrences){
		return DigestSuffix.bytes(hbmFamilyListBytes(repIds,occurrences));
	}

	/** Returns the exact UTF-8 bytes whose digest is stored in HBM provenance row 10. */
	static byte[] hbmFamilyListBytes(final String[] repIds, final int[] occurrences){
		if(repIds==null || occurrences==null || repIds.length!=occurrences.length){
			throw new IllegalArgumentException("HBM family-list arrays must have equal non-null lengths.");
		}
		final StringBuilder text=new StringBuilder(32+repIds.length*48);
		text.append("#rank\trep_id\tocc_total\n");
		for(int rank=0; rank<repIds.length; rank++){
			if(repIds[rank]==null || repIds[rank].length()==0 || occurrences[rank]<1){
				throw new IllegalArgumentException("Invalid HBM family-list projection at rank "+rank+".");
			}
			text.append(rank).append('\t').append(repIds[rank]).append('\t')
				.append(occurrences[rank]).append('\n');
		}
		return text.toString().getBytes(StandardCharsets.UTF_8);
	}

	/** Reads the metadata and 33-column rows without interpreting thresholds. */
	private static Parsed parse(final String path){
		final HashMap<String,String> header=new HashMap<String,String>();
		final ArrayList<String[]> rows=new ArrayList<String[]>();
		final ByteFile input=ByteFile.makeByteFile(path,false);
		boolean columns=false;
		try{
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0){throw new IllegalArgumentException("Blank schema-7 threshold line: "+path);}
				final String text=new String(line,StandardCharsets.US_ASCII);
				if(line[0]=='#'){
					if(columns){throw new IllegalArgumentException("Late schema-7 threshold metadata: "+path);}
					final int tab=text.indexOf('\t');
					if(tab<2 || tab==text.length()-1 || header.put(text.substring(1,tab),text.substring(tab+1))!=null){
						throw new IllegalArgumentException("Malformed or duplicate schema-7 threshold metadata: "+text);
					}
				}else if(!columns){
					columns=true;
				}else{
					final String[] fields=text.split("\\t",-1);
					if(fields.length!=33){throw new IllegalArgumentException("Schema-7 threshold row width "+fields.length+" != 33 in "+path);}
					rows.add(fields);
				}
			}
		}finally{input.close();}
		if(!columns){throw new IllegalArgumentException("Schema-7 threshold table is missing: "+path);}
		return new Parsed(header,rows);
	}

	/** Runs the maintained strict verifier against the artifact's sealed metadata. */
	private static void verifyArtifact(final String path, final String sha80,
			final HashMap<String,String> header){
		final Schema7ThresholdVerifier.Config c=new Schema7ThresholdVerifier.Config();
		c.in=path; c.inSha80=sha80;
		c.expectedFamilies=positiveInt(header,"families",path);
		c.expectedMembers=positiveLong(header,"members",path);
		c.tableSha80=hash(header,"threshold_table_sha80",path);
		c.rosterSha80=hash(header,"roster_sha80",path);
		c.baseLandmarksSha80=hash(header,"base_landmarks_sha80",path);
		c.baseProvenanceSha80=hash(header,"base_provenance_sha80",path);
		c.lengthBoundsSha80=hash(header,"length_bounds_sha80",path);
		c.finalEmpiricalSha80=hash(header,"final_empirical_sha80",path);
		c.consensusSha80=hash(header,"consensus_sha80",path);
		c.coreSha80=hash(header,"core_sha80",path);
		c.coveringSetsSha80=hash(header,"covering_sets_sha80",path);
		c.sidecarSha80=hash(header,"sidecar_sha80",path);
		c.roleManifestSha80=hash(header,"role_manifest_sha80",path);
		c.hbmBundleSha80=hash(header,"hbm_bundle_sha80",path);
		c.hbmProvenanceSha80=hash(header,"hbm_provenance_sha80",path);
		c.d125OutputSha80=hash(header,"d125_output_manifest_sha80",path);
		c.d125CurveSha80=hash(header,"d125_hbm_path_curve_sha80",path);
		c.d125SummarySha80=hash(header,"d125_floor_summary_sha80",path);
		c.d127IdentityCurveSha80=hash(header,"d127_identity_curve_sha80",path);
		c.d127RCurveSha80=hash(header,"d127_r_curve_sha80",path);
		c.d127ComparisonSha80=hash(header,"d127_comparison_manifest_sha80",path);
		Schema7ThresholdVerifier.run(c);
	}

	private static void checkHash(final String path, final String expected,
			final String label){
		final String observed=MagQCTextResource.sha80(path);
		if(!expected.equals(observed)){
			throw new IllegalArgumentException("Schema-7 "+label+" content hash "+observed+" != "+expected+": "+path);
		}
	}

	/** HBM pins identify stored transport bytes; table pins identify decoded text. */
	private static void checkStoredHash(final String path, final String expected, final String label){
		final String observed=DigestSuffix.file(path);
		if(!expected.equals(observed)){
			throw new IllegalArgumentException("Schema-7 "+label+" hash "+observed+
				" != "+expected+": "+path);
		}
	}

	private static String required(final HashMap<String,String> header,
			final String key, final String path){
		final String value=header.get(key);
		if(value==null || value.length()==0){throw new IllegalArgumentException("Missing #"+key+" in "+path);}
		return value;
	}

	private static String hash(final HashMap<String,String> header,
			final String key, final String path){
		return DigestSuffix.requireSuffix(required(header,key,path),key);
	}

	private static int positiveInt(final HashMap<String,String> header,
			final String key, final String path){
		final long value=positiveLong(header,key,path);
		if(value>Integer.MAX_VALUE){throw new IllegalArgumentException("#"+key+" exceeds int range in "+path);}
		return (int)value;
	}

	private static long positiveLong(final HashMap<String,String> header,
			final String key, final String path){
		final String text=required(header,key,path);
		final long value;
		try{value=Long.parseLong(text);}catch(NumberFormatException e){throw new IllegalArgumentException("Non-integer #"+key+" in "+path,e);}
		if(value<1 || !text.equals(Long.toString(value))){throw new IllegalArgumentException("Noncanonical positive #"+key+" in "+path);}
		return value;
	}

	private static int canonicalInt(final String text, final String label,
			final int rank, final String path){
		final int value;
		try{value=Integer.parseInt(text);}catch(NumberFormatException e){throw new IllegalArgumentException("Non-integer "+label+" at rank "+rank+" in "+path,e);}
		if(!text.equals(Integer.toString(value))){throw new IllegalArgumentException("Noncanonical "+label+" at rank "+rank+" in "+path);}
		return value;
	}

	private static double canonicalDouble(final String text, final String label,
			final int rank, final String path){
		final double value;
		try{value=Double.parseDouble(text);}catch(NumberFormatException e){throw new IllegalArgumentException("Non-numeric "+label+" at rank "+rank+" in "+path,e);}
		if(!Double.isFinite(value) || !text.equals(Double.toString(value))){throw new IllegalArgumentException("Noncanonical "+label+" at rank "+rank+" in "+path);}
		return value;
	}

	private static float canonicalFloat(final String text, final String label,
			final int rank, final String path){
		final float value;
		try{value=Float.parseFloat(text);}catch(NumberFormatException e){throw new IllegalArgumentException("Non-numeric "+label+" at rank "+rank+" in "+path,e);}
		if(!Float.isFinite(value) || !text.equals(Float.toString(value))){throw new IllegalArgumentException("Noncanonical "+label+" at rank "+rank+" in "+path);}
		return value;
	}

	private static boolean canonicalBoolean(final String text, final String label,
			final int rank, final String path){
		if("true".equals(text)){return true;}
		if("false".equals(text)){return false;}
		throw new IllegalArgumentException("Noncanonical "+label+" at rank "+rank+" in "+path);
	}

	private static final class Parsed {
		final HashMap<String,String> header;
		final List<String[]> rows;
		Parsed(final HashMap<String,String> header_, final List<String[]> rows_){header=header_; rows=rows_;}
	}

	private static final class RuntimeRoster {
		final int[] familyIds,members;
		final String[] repIds,status;
		final String coreFamilyListSha80;
		RuntimeRoster(final int[] familyIds_, final String[] repIds_,
				final int[] members_, final String[] status_,
				final String coreFamilyListSha80_){
			familyIds=familyIds_; repIds=repIds_; members=members_; status=status_;
			coreFamilyListSha80=coreFamilyListSha80_;
		}
	}

	/** Complete verified schema-7 table and its live-resource bindings. */
	static final class Loaded {
		final String artifactSha80,tableSha80,rosterSha80,coreFamilyListSha80;
		final String consensusSha80;
		final String roleManifestSha80,coreSha80,coveringSetsSha80,sidecarSha80;
		final String hbmBundleSha80,hbmProvenanceSha80,coveringAlphabet;
		final int coveringK;
		final long members;
		final int[] familyIds,n,lengthLo,lengthHi,minRaw,minKmer;
		final String[] repIds,calibrationStatus;
		final boolean[] enabled;
		final double[] minIdentity,minR,minMutual,minKmerDensity;
		final float[] minHbmPath;
		final FamilyRoleManifest roles;
		Loaded(final String artifactSha80_, final String tableSha80_,
				final String rosterSha80_, final String coreFamilyListSha80_,
				final String consensusSha80_,
				final String roleManifestSha80_, final String coreSha80_,
				final String coveringSetsSha80_, final String sidecarSha80_,
				final String hbmBundleSha80_, final String hbmProvenanceSha80_,
				final String coveringAlphabet_, final int coveringK_,
				final long members_, final int[] familyIds_, final String[] repIds_,
				final int[] n_, final String[] calibrationStatus_,
				final boolean[] enabled_, final int[] lengthLo_,
				final int[] lengthHi_, final int[] minRaw_,
				final double[] minIdentity_, final double[] minR_,
				final double[] minMutual_, final int[] minKmer_,
				final double[] minKmerDensity_, final float[] minHbmPath_,
				final FamilyRoleManifest roles_){
			artifactSha80=artifactSha80_; tableSha80=tableSha80_;
			rosterSha80=rosterSha80_; coreFamilyListSha80=coreFamilyListSha80_;
			consensusSha80=consensusSha80_;
			roleManifestSha80=roleManifestSha80_; coreSha80=coreSha80_;
			coveringSetsSha80=coveringSetsSha80_; sidecarSha80=sidecarSha80_;
			hbmBundleSha80=hbmBundleSha80_; hbmProvenanceSha80=hbmProvenanceSha80_;
			coveringAlphabet=coveringAlphabet_; coveringK=coveringK_; members=members_;
			familyIds=familyIds_; repIds=repIds_; n=n_;
			calibrationStatus=calibrationStatus_; enabled=enabled_;
			lengthLo=lengthLo_; lengthHi=lengthHi_; minRaw=minRaw_;
			minIdentity=minIdentity_; minR=minR_; minMutual=minMutual_;
			minKmer=minKmer_; minKmerDensity=minKmerDensity_;
			minHbmPath=minHbmPath_; roles=roles_;
		}
	}

	private static final String ROSTER_COLUMNS=
		"active_index\tfamily_id\trep_id\tlineage_provenance\tparent_family_id"+
		"\tchild_index\tmember_count\tsequence_bytes\tmember_file\tmember_sha80"+
		"\tsource_member_path\tsource_member_sha80\tsource_consensus_path"+
		"\tconsensus_sequence_sha80\tsource_consensus_artifact_sha80"+
		"\tmember_provenance\tretained_n\tcalibration_status";
	private static final String ROSTER_SCHEMA="post_split_final_roster_v4";
}
