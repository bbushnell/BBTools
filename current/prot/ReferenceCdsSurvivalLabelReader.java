package prot;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.HashSet;

/**
 * Strict consumer for one selection-specific ReferenceCdsSurvivalSidecar.
 * This reader computes labels from original reference identities only; cache
 * family counts, fresh calls, and RNA observations are deliberately invisible.
 */
public final class ReferenceCdsSurvivalLabelReader {
	/** Accepts v1 (raw containing ids) and v2 (IdentifierCodec-encoded containing ids, decoded strictly here). */
	private static final String SCHEMA_V1=ReferenceCdsSurvivalSidecar.SCHEMA_V1, SCHEMA_V2=ReferenceCdsSurvivalSidecar.SCHEMA;
	private static final String HEADER="role\tassembly_id\tgene_key\tsource_contig\tgene_id\tstart1\tend1\tstrand\tstatus\tcontaining_shred_count\tcontaining_shred_ids";
	private ReferenceCdsSurvivalLabelReader(){ }

	/** Four count operands plus the two D39/METHODS gene-survival labels. */
	public static final class Counts {
		public final int nativeTotal, nativeRetained, foreignRetained, binTotalGenes;
		public final double completeness, contamination;
		Counts(int nativeTotal_, int nativeRetained_, int foreignRetained_){
			nativeTotal=nativeTotal_; nativeRetained=nativeRetained_; foreignRetained=foreignRetained_;
			binTotalGenes=nativeRetained+foreignRetained;
			completeness=nativeTotal==0 ? 0d : nativeRetained/(double)nativeTotal;
			contamination=binTotalGenes==0 ? 0d : foreignRetained/(double)binTotalGenes;
		}
	}

	/** Loads, validates, and scores a sidecar.  Expected hashes may be null or "-" for a
	 * caller that has already validated the binding in its own index. */
	public static Counts load(String sidecarPath, String expectedSidecarSha256,
			String expectedSourceManifestSha256) {
		if(sidecarPath==null || sidecarPath.length()==0){throw new IllegalArgumentException("Missing sidecar path");}
		if(expectedSidecarSha256!=null && !expectedSidecarSha256.equals("-")){
			String actual=sha256File(sidecarPath);
			if(!expectedSidecarSha256.equalsIgnoreCase(actual)){
				throw new IllegalArgumentException("Sidecar SHA-256 mismatch for "+sidecarPath+": expected "+expectedSidecarSha256+" got "+actual);
			}
		}
		try{
			return read(sidecarPath, expectedSourceManifestSha256);
		}catch(IOException e){throw new RuntimeException("Could not read reference-CDS sidecar "+sidecarPath,e);}
	}

	private static Counts read(String file, String expectedSourceManifestSha256) throws IOException{
		boolean schema=false, sourceManifest=false, header=false, encoded=false;
		String sourceHash=null;
		int nativeTotal=0, nativeRetained=0, foreignRetained=0;
		HashSet<String> keys=new HashSet<String>();
		try(BufferedReader r=Files.newBufferedReader(Paths.get(file), StandardCharsets.UTF_8)){
			String line; int lineNo=0;
			while((line=r.readLine())!=null){
				lineNo++;
				if(line.isEmpty()){continue;}
				if(line.charAt(0)=='#'){
					if(line.startsWith("#schema_version\t")){
						if(schema){throw bad(file,lineNo,"duplicate schema");}
						if(line.equals("#schema_version\t"+SCHEMA_V2)){encoded=true;}
						else if(!line.equals("#schema_version\t"+SCHEMA_V1)){throw bad(file,lineNo,"invalid schema");}
						schema=true;
					}else if(line.startsWith("#source_manifest_sha256\t")){
						if(sourceManifest){throw bad(file,lineNo,"duplicate source-manifest hash");}
						sourceHash=line.substring("#source_manifest_sha256\t".length());
						if(sourceHash.length()==0){throw bad(file,lineNo,"empty source-manifest hash");}
						sourceManifest=true;
					}
					continue;
				}
				if(!header){
					if(!line.equals(HEADER)){throw bad(file,lineNo,"invalid data header");}
					header=true; continue;
				}
				String[] f=line.split("\t",-1);
				if(f.length!=11){throw bad(file,lineNo,"expected 11 fields, got "+f.length);}
				if(!(f[0].equals("native") || f[0].equals("foreign"))){throw bad(file,lineNo,"invalid role "+f[0]);}
				if(f[1].isEmpty() || f[2].isEmpty() || f[3].isEmpty() || f[4].isEmpty()){throw bad(file,lineNo,"empty identity field");}
				if(!keys.add(f[2])){throw bad(file,lineNo,"duplicate gene_key "+f[2]);}
				final String status=f[8];
				if(!(status.equals("RETAINED") || status.equals("SPLIT_OR_PARTIAL") || status.equals("OMITTED"))){throw bad(file,lineNo,"invalid status "+status);}
				final int containing=parseNonnegative(f[9],file,lineNo,"containing_shred_count");
				if(containing==0){if(!f[10].equals("-")){throw bad(file,lineNo,"zero containing count with nonempty IDs");}}
				else{
					if(f[10].equals("-")){throw bad(file,lineNo,"positive containing count with empty IDs");}
					String[] ids=f[10].split(",",-1); if(ids.length!=containing){throw bad(file,lineNo,"containing count does not match IDs");}
					HashSet<String> seen=new HashSet<String>();
					for(String id0 : ids){
						String id=id0;
						if(encoded){try{id=IdentifierCodec.decode(id0);}catch(IllegalArgumentException e){throw bad(file,lineNo,"malformed encoded containing shred ID: "+e.getMessage());}}
						if(id.isEmpty() || !seen.add(id)){throw bad(file,lineNo,"duplicate/empty containing shred ID");}
					}
				}
				if(f[0].equals("native")){nativeTotal++; if(status.equals("RETAINED")){nativeRetained++;}}
				else if(status.equals("RETAINED")){foreignRetained++;}
			}
		}
		if(!schema || !sourceManifest || !header){throw new IllegalArgumentException("Incomplete reference-CDS sidecar headers: "+file);}
		if(expectedSourceManifestSha256!=null && !expectedSourceManifestSha256.equals("-")
			&& !expectedSourceManifestSha256.equalsIgnoreCase(sourceHash)){
			throw new IllegalArgumentException("Source-manifest SHA-256 mismatch for "+file+": expected "+expectedSourceManifestSha256+" got "+sourceHash);
		}
		return new Counts(nativeTotal,nativeRetained,foreignRetained);
	}

	private static int parseNonnegative(String s,String file,int line,String field){
		try{int v=Integer.parseInt(s); if(v<0){throw new NumberFormatException();} return v;}
		catch(NumberFormatException e){throw bad(file,line,"invalid "+field+": "+s);}
	}
	private static IllegalArgumentException bad(String file,int line,String msg){return new IllegalArgumentException(file+":"+line+": "+msg);}
	public static String sha256File(String file){
		try{
			MessageDigest md=MessageDigest.getInstance("SHA-256"); byte[] buf=new byte[1<<16];
			try(java.io.InputStream in=Files.newInputStream(Paths.get(file))){for(int n; (n=in.read(buf))>=0;){if(n>0){md.update(buf,0,n);}}}
			StringBuilder out=new StringBuilder(64); for(byte b : md.digest()){out.append(String.format("%02x", b));} return out.toString();
		}catch(Exception e){throw new RuntimeException("Could not hash "+file,e);}
	}
}
