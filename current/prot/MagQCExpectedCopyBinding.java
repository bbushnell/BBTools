package prot;

import java.nio.charset.StandardCharsets;
import java.util.LinkedHashMap;
import fileIO.ByteFile;

/**
 * Binds one hash-verified {@link MagQCExpectedCopyTableTyped.Table} to one hash-verified
 * {@link MagQCNetBundle} through a synthetic declaration file that names, for each bundle
 * subnet, which table population ({@code global} or one {@code phylum}) supplies its frozen
 * expectations. The binding precompiles one {@link MagQCExpectedCopyFeatures.Mapping} per
 * bundle subnet, IN BUNDLE ORDER, once at {@link #bind} time; every later {@link #compute}
 * call only replays that frozen mapping.
 *
 * <p>The binding never consults a bin's own phylum and never derives anything from a
 * subnet's name: the population a subnet draws from comes exclusively from the declaration
 * row keyed by the subnet's id. Two subnets with different names can be bound to the same
 * population, and one subnet's declared population can be swapped (declaration A vs B)
 * without touching the bundle, the table, or the subnet's name.
 *
 * <p>Retains four validated identities for downstream provenance: the verified table hash,
 * the verified declaration hash, the family-list hash shared by the table and the bundle,
 * and the release-manifest hash shared by the bundle and the declaration.
 */
public final class MagQCExpectedCopyBinding {

	private static final String DECLARATION_SCHEMA="subnet_population_declaration_v1";
	private static final String DECLARATION_COLUMNS="subnet_id\tpopulation_type\tpopulation_name";

	private final String tableSha256,declarationSha256,familyListSha256,releaseManifestSha256;
	private final int numFam, itemWidth;
	private final String[] subnetIds,populationTypes,populationNames;
	private final MagQCExpectedCopyFeatures.Mapping[] mappings;

	private MagQCExpectedCopyBinding(String tableSha256,String declarationSha256,String familyListSha256,
			String releaseManifestSha256,int numFam,int itemWidth,String[] subnetIds,String[] populationTypes,
			String[] populationNames,MagQCExpectedCopyFeatures.Mapping[] mappings){
		this.tableSha256=tableSha256;
		this.declarationSha256=declarationSha256;
		this.familyListSha256=familyListSha256;
		this.releaseManifestSha256=releaseManifestSha256;
		this.numFam=numFam;
		this.itemWidth=itemWidth;
		this.subnetIds=subnetIds.clone();
		this.populationTypes=populationTypes.clone();
		this.populationNames=populationNames.clone();
		this.mappings=mappings.clone();
	}

	/**
	 * Verifies the table hash, the family-list triangle, the canonical item order, the
	 * declaration hash and its release-manifest/id-set agreement with the bundle, then
	 * compiles one {@link MagQCExpectedCopyFeatures.Mapping} per bundle subnet in bundle
	 * order. Fails closed with {@link IllegalArgumentException} on the first violation.
	 */
	public static MagQCExpectedCopyBinding bind(String tablePath,String expectedTableSha256,MagQCNetBundle bundle,
			String declarationPath,String expectedDeclarationSha256,int numFam){
		if(bundle==null){fail("bundle is required");}
		if(numFam<=0){fail("numFam must be positive");}
		final MagQCExpectedCopyTableTyped.Table table=MagQCExpectedCopyTableTyped.open(tablePath,expectedTableSha256);

		final String fl=bundle.metadata("family_list_sha256");
		if(!matchesBundleDigest(bundle,fl,table.familyListHash)){fail("family list hash mismatch between table and bundle");}

		final int width=table.items.length, nonprotein=width-numFam;
		if(nonprotein!=MagQCObservationLayout.LEGACY_COUNT && nonprotein!=MagQCObservationLayout.EXTENDED_COUNT){
			fail("table item width must match numFam+5 or numFam+70");
		}
		for(int i=0;i<numFam;i++){
			final MagQCExpectedCopyTableTyped.Item it=table.items[i];
			if(!"P".equals(it.type)||!it.key.equals(Integer.toString(i))||it.ordinal!=i){
				fail("canonical item order mismatch at "+i);
			}
		}
		for(int i=numFam;i<table.items.length;i++){
			final MagQCExpectedCopyTableTyped.Item it=table.items[i];
			final String expectedKey=MagQCObservationLayout.nonproteinKey(i-numFam);
			if(!"N".equals(it.type)||!it.key.equals(expectedKey)||it.ordinal!=i){
				fail("canonical item order mismatch at "+i);
			}
		}

		if(!MagQCExpectedCopyTableTyped.validDigest(expectedDeclarationSha256)){
			fail("invalid expected declaration hash");
		}
		final String actualDeclarationSha256=HbmMemberIndexFormat.toHexLower(MagQCTextResource.digest(declarationPath));
		if(!MagQCExpectedCopyTableTyped.matchesDigest(expectedDeclarationSha256,actualDeclarationSha256)){
			fail("declaration hash mismatch");
		}
		final Declaration declaration=parseDeclaration(declarationPath);

		final String rm=bundle.metadata("release_manifest_sha256");
		if(!matchesBundleDigest(bundle,rm,declaration.releaseManifestSha256)){
			fail("release manifest hash mismatch between bundle and declaration");
		}

		for(int i=0;i<bundle.size();i++){
			final String id=bundle.subnet(i).id;
			if(!declaration.rows.containsKey(id)){fail("declaration missing subnet id "+id);}
		}
		if(declaration.rows.size()!=bundle.size()){
			for(String id:declaration.rows.keySet()){
				boolean found=false;
				for(int i=0;i<bundle.size();i++){if(bundle.subnet(i).id.equals(id)){found=true;break;}}
				if(!found){fail("declaration has extra subnet id "+id);}
			}
			fail("declaration/bundle subnet id set mismatch");
		}

		final String[] subnetIds=new String[bundle.size()];
		final String[] populationTypes=new String[bundle.size()];
		final String[] populationNames=new String[bundle.size()];
		final MagQCExpectedCopyFeatures.Mapping[] mappings=new MagQCExpectedCopyFeatures.Mapping[bundle.size()];
		for(int i=0;i<bundle.size();i++){
			final MagQCNetBundle.Subnet s=bundle.subnet(i);
			final DeclarationRow row=declaration.rows.get(s.id);
			final MagQCExpectedCopyTableTyped.Population pop=table.population(row.type,row.name);
			final int[] indexes;
			final int specialCount=MagQCObservationLayout.specialCount(s.type);
			if(specialCount>=0){
				if(s.familyRanks!=null || s.numObs!=specialCount){
					fail("special subnet observation layout mismatch: "+s.id+" type="+s.type);
				}
				final int offset="trna_anticodon".equals(s.type) ? MagQCObservationLayout.LEGACY_COUNT : 0;
				if(offset+specialCount>nonprotein){fail("trna_anticodon subnet requires an extended expected-copy table: "+s.id);}
				indexes=new int[specialCount];
				for(int j=0; j<indexes.length; j++){indexes[j]=numFam+offset+j;}
			}else{
				if(s.familyRanks==null){fail("famset subnet requires family ranks: "+s.id);}
				if(s.familyRanks.length!=s.numObs){fail("family rank count does not match numObs: "+s.id);}
				for(int r:s.familyRanks){
					if(r<0||r>=numFam){fail("family rank out of range: "+s.id+" rank "+r);}
				}
				indexes=s.familyRanks.clone();
			}
			mappings[i]=MagQCExpectedCopyFeatures.compileMapping(width,indexes,pop.means());
			subnetIds[i]=s.id;
			populationTypes[i]=row.type;
			populationNames[i]=row.name;
		}

		return new MagQCExpectedCopyBinding(expectedTableSha256,expectedDeclarationSha256,fl,rm,numFam,width,
				subnetIds,populationTypes,populationNames,mappings);
	}

	/**
	 * Compares a verified legacy artifact's full digest using the bundle's declared algorithm.
	 * The table/declaration byte checks retain all their original digest bits. Only this
	 * cross-format comparison projects a full digest to sha80, and only for a sha80 bundle.
	 */
	private static boolean matchesBundleDigest(MagQCNetBundle bundle,String stored,String legacy){
		assert(bundle!=null) : "Digest comparison requires the validated bundle's hash algorithm";
		if(!MagQCExpectedCopyTableTyped.validDigest(legacy)){fail("invalid artifact digest");}
		final String algorithm=bundle.metadata("hash_algorithm");
		if("SHA-256".equals(algorithm)){
			if(stored==null || !stored.matches("[0-9a-f]{64}")){fail("invalid full bundle digest");}
			return stored.equals(legacy);
		}
		if("SHA-256/80".equals(algorithm)){
			if(stored==null || !stored.matches("[0-9a-f]{20}")){fail("invalid sha80 bundle digest");}
			return stored.equals(legacy.substring(legacy.length()-20));
		}
		fail("unsupported bundle hash algorithm");
		return false;
	}

	public int size(){return mappings.length;}
	public int numFam(){return numFam;}
	/** Returns the bound table's actual width, including its legacy or extended N block. */
	public int itemWidth(){return itemWidth;}
	public String subnetId(int order){checkOrder(order);return subnetIds[order];}
	public String populationType(int order){checkOrder(order);return populationTypes[order];}
	public String populationName(int order){checkOrder(order);return populationNames[order];}
	public String tableSha256(){return tableSha256;}
	public String declarationSha256(){return declarationSha256;}
	public String familyListSha256(){return familyListSha256;}
	public String releaseManifestSha256(){return releaseManifestSha256;}

	/** Computes this subnet's frozen features from whole-bin observations. */
	public MagQCExpectedCopyFeatures.Result compute(int order,int[] observedByItem){
		checkOrder(order);
		return mappings[order].compute(observedByItem);
	}

	private void checkOrder(int order){if(order<0||order>=mappings.length){fail("subnet order out of range: "+order);}}

	private static final class DeclarationRow{final String type,name;DeclarationRow(String t,String n){type=t;name=n;}}
	private static final class Declaration{final String releaseManifestSha256;final LinkedHashMap<String,DeclarationRow> rows;
		Declaration(String r,LinkedHashMap<String,DeclarationRow> w){releaseManifestSha256=r;rows=w;}}

	/** Parses the synthetic {@code subnet_population_declaration_v1} format (see the design brief §2). */
	private static Declaration parseDeclaration(String file){
		String schema=null,releaseHash=null,columns=null;
		final LinkedHashMap<String,DeclarationRow> rows=new LinkedHashMap<String,DeclarationRow>();
		boolean sawData=false;
		final ByteFile input=ByteFile.makeByteFile(file,false);
		try{
			for(byte[] raw=input.nextLine();raw!=null;raw=input.nextLine()){
				final String line=new String(raw,StandardCharsets.UTF_8);
				if(line.length()==0){continue;}
				if(line.charAt(0)=='#'){
					if(sawData){fail("declaration metadata after data");}
					if(line.startsWith("#schema_version\t")){
						final String[] p=line.split("\t",-1);
						if(p.length!=2||schema!=null||!DECLARATION_SCHEMA.equals(p[1])){fail("invalid declaration schema");}
						schema=p[1];
					}else if(line.startsWith("#release_manifest_sha256\t")){
						final String[] p=line.split("\t",-1);
						if(p.length!=2||releaseHash!=null||!MagQCExpectedCopyTableTyped.validDigest(p[1])){fail("invalid declaration release manifest metadata");}
						releaseHash=p[1];
					}else if(line.startsWith("#columns\t")){
						final String c=line.substring(9);
						if(columns!=null||!DECLARATION_COLUMNS.equals(c)){fail("invalid declaration columns");}
						columns=c;
					}else{
						fail("unknown declaration metadata");
					}
				}else{
					sawData=true;
					final String[] p=line.split("\t",-1);
					if(p.length!=3){fail("invalid declaration row field count");}
					final String id=p[0],type=p[1],name=p[2];
					if(!id.matches("[A-Za-z0-9_.-]+")){fail("invalid declaration subnet id "+id);}
					if(!"global".equals(type)&&!"phylum".equals(type)){fail("invalid declaration population type "+type);}
					if("global".equals(type)){
						if(!"-".equals(name)){fail("declaration global row requires name '-': "+id);}
					}else{
						if(name.length()==0||"-".equals(name)||name.indexOf('\t')>=0||name.indexOf('\r')>=0||name.indexOf('\n')>=0){
							fail("invalid declaration phylum name for "+id);
						}
					}
					if(rows.put(id,new DeclarationRow(type,name))!=null){fail("duplicate declaration subnet id "+id);}
				}
			}
		}finally{if(input.close()){fail("I/O error reading population declarations: "+file);}}
		if(schema==null||releaseHash==null||columns==null){fail("incomplete declaration headers");}
		if(rows.isEmpty()){fail("empty declaration");}
		return new Declaration(releaseHash,rows);
	}

	private static void fail(String message){throw new IllegalArgumentException(message);}
}
