package prot;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Consumer-side reader for Sayu's selection-specific reference-CDS index.
 * The index key is (bin-manifest SHA-256, bin_id); no cache aggregate fallback
 * is permitted.  This class is intentionally independent of the producer so
 * a malformed or stale index fails before vector emission.
 */
public final class MagQCReferenceCdsIndex {
	private static final String SCHEMA_V1=ReferenceCdsSidecarIndex.SCHEMA_V1;
	private static final String SCHEMA=ReferenceCdsSidecarIndex.SCHEMA;
	private static final String COLUMNS="bin_manifest_sha256\tbin_id\tsidecar_path\tsidecar_sha256\tsource_manifest_path\tsource_manifest_sha256\tselection_binding_sha256\tnative_source_ids\tforeign_source_ids\tselected_source_span_ids";
	private static final String HEX="[0-9a-fA-F]{64}";
	private final String binManifestPath, binManifestSha256;
	private final int schemaVersion;
	private final HashMap<String,Binding> rows;
	private MagQCReferenceCdsIndex(String path,String hash,int schemaVersion_,HashMap<String,Binding> rows_){
		binManifestPath=path;binManifestSha256=hash;schemaVersion=schemaVersion_;rows=rows_;
	}

	public static final class Binding {
		public final String binManifestSha256, binId, sidecarPath, sidecarSha256,
			sourceManifestPath, sourceManifestSha256, selectionBindingSha256,
			nativeSourceIds, foreignSourceIds, selectedSourceSpanIds;
		Binding(String bm,String id,String sp,String sh,String sm,String smh,String sb,
				String n,String f,String spans){binManifestSha256=bm;binId=id;sidecarPath=sp;sidecarSha256=sh;
			sourceManifestPath=sm;sourceManifestSha256=smh;selectionBindingSha256=sb;nativeSourceIds=n;
			foreignSourceIds=f;selectedSourceSpanIds=spans;}
	}

	public static MagQCReferenceCdsIndex load(String indexPath, String binManifestPath){
		try{
			String actualManifestHash=sha256File(binManifestPath);
			HashMap<String,Binding> rows=new HashMap<String,Binding>();
			String declaredPath=null, declaredHash=null, schemaName=null; boolean schema=false, columns=false;
			try(BufferedReader r=Files.newBufferedReader(Paths.get(indexPath), StandardCharsets.UTF_8)){
				String line; int lineNo=0;
				while((line=r.readLine())!=null){lineNo++; if(line.isEmpty()){continue;}
					if(line.charAt(0)=='#'){
						if(line.startsWith("#schema_version\t")){
							if(schema){throw bad(indexPath,lineNo,"duplicate schema");}
							schemaName=line.substring("#schema_version\t".length());
							if(!(SCHEMA_V1.equals(schemaName) || SCHEMA.equals(schemaName))){throw bad(indexPath,lineNo,"invalid schema");}
							schema=true;
						}
						else if(line.startsWith("#bin_manifest\t")){if(declaredPath!=null){throw bad(indexPath,lineNo,"duplicate bin manifest header");} String[] f=line.split("\t",-1); if(f.length!=5 || !f[2].equals("#") || !f[3].equals("sha256") || !f[4].matches(HEX)){throw bad(indexPath,lineNo,"invalid bin manifest header");} declaredPath=f[1]; declaredHash=f[4];}
						else if(line.startsWith("#columns\t")){if(columns || !line.equals("#columns\t"+COLUMNS)){throw bad(indexPath,lineNo,"invalid columns header");} columns=true;}
						continue;
					}
					if(!columns){throw bad(indexPath,lineNo,"data precedes columns header");}
						String[] f=line.split("\t",-1); if(f.length!=10){throw bad(indexPath,lineNo,"expected 10 fields, got "+f.length);}
						if(!f[0].matches(HEX) || f[1].isEmpty() || f[2].isEmpty() || !f[3].matches(HEX) || f[4].isEmpty() || !f[5].matches(HEX) || !f[6].matches(HEX)){throw bad(indexPath,lineNo,"invalid key/path/hash field");}
						if(!f[0].equalsIgnoreCase(actualManifestHash)){throw bad(indexPath,lineNo,"row bin-manifest hash differs from the bound manifest");}
					for(int i=7;i<10;i++){if(f[i].isEmpty() || f[i].indexOf(',')>=0 || f[i].indexOf('\n')>=0 || f[i].indexOf('\r')>=0){throw bad(indexPath,lineNo,"invalid source binding field");}}
					Binding b=new Binding(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9]);
					if(rows.put(key(f[0],f[1]),b)!=null){throw bad(indexPath,lineNo,"duplicate key");}
				}
			}
			if(!schema || !columns || declaredPath==null){throw new IllegalArgumentException("Incomplete reference-CDS index headers: "+indexPath);}
			if(!declaredPath.equals(binManifestPath) || !declaredHash.equalsIgnoreCase(actualManifestHash)){throw new IllegalArgumentException("Bin-manifest binding mismatch: index="+declaredPath+"/"+declaredHash+" actual="+binManifestPath+"/"+actualManifestHash);}
			final MagQCBinManifest.Manifest manifest=MagQCBinManifest.load(binManifestPath);
			final int indexSchema=SCHEMA.equals(schemaName) ? MagQCBinManifest.SCHEMA_VERSION_ENCODED
				: MagQCBinManifest.SCHEMA_VERSION;
			if(manifest.schemaVersion!=indexSchema){throw new IllegalArgumentException("Index schema does not match bound bin manifest: index="+schemaName+" manifest="+manifest.schemaVersion);}
			final HashSet<String> manifestBinIds=new HashSet<String>();
			for(MagQCBinManifest.Bin b : manifest.bins){manifestBinIds.add(b.id);}
			for(Binding b : rows.values()){
				if(!manifestBinIds.contains(b.binId)){throw new IllegalArgumentException("Index row bin_id absent from bound manifest: "+b.binId);}
			}
			return new MagQCReferenceCdsIndex(binManifestPath,actualManifestHash,indexSchema,rows);
		}catch(IOException e){throw new RuntimeException("Could not read reference-CDS index "+indexPath,e);}
	}

	/** Resolves and validates one row against the exact recorded bin selection. */
	public ReferenceCdsSurvivalLabelReader.Counts labelsFor(MagQCBinManifest.Bin bin){
		if(bin.schemaVersion!=schemaVersion){throw new IllegalArgumentException("Reference-CDS bin schema mismatch for "+bin.id);}
		Binding b=rows.get(key(binManifestSha256,bin.id));
		if(b==null){throw new IllegalArgumentException("No reference-CDS index row for key "+binManifestSha256+","+bin.id);}
		if(!bin.nativeContigs.equals(b.nativeSourceIds) || !bin.foreignContigs.equals(b.foreignSourceIds)){
			throw new IllegalArgumentException("Reference-CDS selection binding mismatch for "+bin.id);
		}
		// Use the producer's public canonical helper. Span IDs remain a separate
		// explicit index field and are intentionally not part of this hash.
		String actualBinding=ReferenceCdsSidecarIndex.selectionBindingSha256(
			b.nativeSourceIds,b.foreignSourceIds);
		if(!b.selectionBindingSha256.equalsIgnoreCase(actualBinding)){throw new IllegalArgumentException("Selection binding hash mismatch for "+bin.id);}
		if(!b.sourceManifestSha256.equalsIgnoreCase(sha256File(b.sourceManifestPath))){throw new IllegalArgumentException("Source-manifest hash mismatch for "+bin.id);}
		ReferenceCdsSidecarIndex.validateSelectedSourceSpanIds(b.sourceManifestPath, b.sidecarPath,
				bin, b.selectedSourceSpanIds);
		return ReferenceCdsSurvivalLabelReader.load(b.sidecarPath,b.sidecarSha256,b.sourceManifestSha256);
	}

	private static String key(String hash,String id){return hash.toLowerCase()+"\u0000"+id;}
	private static IllegalArgumentException bad(String f,int n,String s){return new IllegalArgumentException(f+":"+n+": "+s);}
	private static String sha256File(String file){try{MessageDigest md=MessageDigest.getInstance("SHA-256"); byte[] buf=new byte[1<<16]; try(java.io.InputStream in=Files.newInputStream(Paths.get(file))){for(int n;(n=in.read(buf))>=0;){if(n>0){md.update(buf,0,n);}}} StringBuilder b=new StringBuilder(64); for(byte x:md.digest()){b.append(String.format("%02x",x));} return b.toString();}catch(Exception e){throw new RuntimeException("Could not hash "+file,e);}}
}
