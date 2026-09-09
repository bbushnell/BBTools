package prot;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.LinkedHashMap;

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
	private static final String HEX="0123456789abcdef";

	private final String tableSha256,declarationSha256,familyListSha256,releaseManifestSha256;
	private final int numFam;
	private final String[] subnetIds,populationTypes,populationNames;
	private final MagQCExpectedCopyFeatures.Mapping[] mappings;

	private MagQCExpectedCopyBinding(String tableSha256,String declarationSha256,String familyListSha256,
			String releaseManifestSha256,int numFam,String[] subnetIds,String[] populationTypes,
			String[] populationNames,MagQCExpectedCopyFeatures.Mapping[] mappings){
		this.tableSha256=tableSha256;
		this.declarationSha256=declarationSha256;
		this.familyListSha256=familyListSha256;
		this.releaseManifestSha256=releaseManifestSha256;
		this.numFam=numFam;
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
		if(fl==null||!fl.equals(table.familyListHash)){fail("family list hash mismatch between table and bundle");}

		if(table.items.length!=numFam+5){fail("table item width does not match numFam+5");}
		for(int i=0;i<numFam;i++){
			final MagQCExpectedCopyTableTyped.Item it=table.items[i];
			if(!"P".equals(it.type)||!it.key.equals(Integer.toString(i))||it.ordinal!=i){
				fail("canonical item order mismatch at "+i);
			}
		}
		for(int i=numFam;i<table.items.length;i++){
			final MagQCExpectedCopyTableTyped.Item it=table.items[i];
			final String expectedKey=MagQCExpectedCopyTableTyped.NCRNA[i-numFam];
			if(!"N".equals(it.type)||!it.key.equals(expectedKey)||it.ordinal!=i){
				fail("canonical item order mismatch at "+i);
			}
		}

		if(expectedDeclarationSha256==null||!expectedDeclarationSha256.matches("[0-9a-f]{64}")){
			fail("invalid expected declaration hash");
		}
		final String actualDeclarationSha256=sha256(declarationPath);
		if(!actualDeclarationSha256.equals(expectedDeclarationSha256)){
			fail("declaration hash mismatch: expected "+expectedDeclarationSha256+" actual "+actualDeclarationSha256);
		}
		final Declaration declaration=parseDeclaration(declarationPath);

		final String rm=bundle.metadata("release_manifest_sha256");
		if(rm==null||!rm.equals(declaration.releaseManifestSha256)){
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

		final int width=numFam+5;
		final String[] subnetIds=new String[bundle.size()];
		final String[] populationTypes=new String[bundle.size()];
		final String[] populationNames=new String[bundle.size()];
		final MagQCExpectedCopyFeatures.Mapping[] mappings=new MagQCExpectedCopyFeatures.Mapping[bundle.size()];
		for(int i=0;i<bundle.size();i++){
			final MagQCNetBundle.Subnet s=bundle.subnet(i);
			final DeclarationRow row=declaration.rows.get(s.id);
			final MagQCExpectedCopyTableTyped.Population pop=table.population(row.type,row.name);
			final int[] indexes;
			if(s.familyRanks==null){
				if(s.numObs!=5){fail("ncrna subnet must declare 5 observations: "+s.id);}
				indexes=new int[]{numFam,numFam+1,numFam+2,numFam+3,numFam+4};
			}else{
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

		return new MagQCExpectedCopyBinding(expectedTableSha256,actualDeclarationSha256,fl,rm,numFam,
				subnetIds,populationTypes,populationNames,mappings);
	}

	public int size(){return mappings.length;}
	public int numFam(){return numFam;}
	public int itemWidth(){return numFam+5;}
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
		try(BufferedReader br=new BufferedReader(new InputStreamReader(new FileInputStream(file),StandardCharsets.UTF_8))){
			for(String line;(line=br.readLine())!=null;){
				if(line.length()==0){continue;}
				if(line.charAt(0)=='#'){
					if(sawData){fail("declaration metadata after data");}
					if(line.startsWith("#schema_version\t")){
						final String[] p=line.split("\t",-1);
						if(p.length!=2||schema!=null||!DECLARATION_SCHEMA.equals(p[1])){fail("invalid declaration schema");}
						schema=p[1];
					}else if(line.startsWith("#release_manifest_sha256\t")){
						final String[] p=line.split("\t",-1);
						if(p.length!=2||releaseHash!=null||!p[1].matches("[0-9a-f]{64}")){fail("invalid declaration release manifest metadata");}
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
		}catch(IOException e){throw new RuntimeException(e);}
		if(schema==null||releaseHash==null||columns==null){fail("incomplete declaration headers");}
		if(rows.isEmpty()){fail("empty declaration");}
		return new Declaration(releaseHash,rows);
	}

	private static String sha256(String file){
		try{
			final MessageDigest md=MessageDigest.getInstance("SHA-256");
			try(FileInputStream in=new FileInputStream(file)){
				final byte[] buf=new byte[65536];
				for(int n;(n=in.read(buf))!=-1;){md.update(buf,0,n);}
			}
			final byte[] d=md.digest();
			final StringBuilder b=new StringBuilder(64);
			for(byte x:d){b.append(HEX.charAt((x>>>4)&15)).append(HEX.charAt(x&15));}
			return b.toString();
		}catch(IOException e){throw new RuntimeException(e);}
		catch(NoSuchAlgorithmException e){throw new RuntimeException(e);}
	}

	private static void fail(String message){throw new IllegalArgumentException(message);}
}
