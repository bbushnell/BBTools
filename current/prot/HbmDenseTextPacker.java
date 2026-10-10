package prot;

import static prot.HbmDenseTextBundle.digest;
import static prot.HbmDenseTextBundle.update;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;

/**
 * Versioned dense HBM text transport. Families have sha80 checksums, and the
 * root checksum covers all header lines and the ordered family-checksum lines.
 * Every digest consumes the UTF-8 line bytes followed by one LF. CRLF and LF
 * therefore have identical semantics under ByteFile's line normalization.
 *
 * Sixteen ordered provenance pins retain the native runtime/build bindings.
 * The graph contract names MQHBv1's canonical knobs, pad0, 22 residues and
 * weight=count. This transport does not authenticate an untrusted model author.
 *
 * @author Collei
 */
public final class HbmDenseTextPacker {

	static final String FORMAT=HbmDenseTextBundle.FORMAT;
	static final String CONTRACT=HbmDenseTextBundle.CONTRACT;

	/**
	 * Exports a validated native bundle to a fresh, uncompressed .hbmt file.
	 * Compression is independent: bgzip -l9 preserves the normalized text contract.
	 * Exact graph comparison against the trusted native loader is required before
	 * this method returns; failures remove only the newly created output.
	 */
	static void pack(final Path input, final String reference, final String provenance,
			final Path output) throws Exception{
		if(Files.exists(output) || !output.toString().endsWith(".hbmt")){
			throw new IllegalArgumentException("A fresh uncompressed .hbmt output is required");
		}
		final List<ProteinSequence> sequences=ProteinSearch.readFasta(reference);
		final ArrayList<String> roster=new ArrayList<String>(sequences.size());
		final HashMap<String, byte[]> consensus=new HashMap<String, byte[]>();
		for(final ProteinSequence sequence : sequences){
			roster.add(sequence.id);
			if(consensus.put(sequence.id, sequence.enc)!=null){throw new IllegalArgumentException("Duplicate consensus: "+sequence.id);}
		}
		final byte[][] trusted=HbmBundleLoader.loadSemanticProvenance(provenance);
		final String originalPin=DigestSuffix.file(input.toString());
		final HbmBundleLoader.Loaded expected=HbmBundleLoader.load(input, roster, consensus::get, trusted);
		final Path temp=Files.createTempDirectory("hbm-dense-pack-"), dump=temp.resolve("dense.txt");
		boolean created=false, complete=false;
		try{
			HbmTextDump.run(input, null, reference, dump, HbmTextDump.Encoding.DENSE, HbmTextDump.CoordMode.EXPLICIT, 1, false);
			if(!originalPin.equals(DigestSuffix.file(input.toString()))){throw new IOException("Native input changed during export");}
			// CREATE_NEW avoids destroying a destination created concurrently.
			Files.createFile(output); created=true;
			envelope(dump.toString(), output.toString(), trusted);
			HbmDenseTextLoader.loadVerified(output.toString(), roster, consensus::get, trusted, 1).assertStructuralEquivalent(expected);
			complete=true;
		}finally{
			Files.deleteIfExists(dump); Files.deleteIfExists(temp);
			if(created && !complete){Files.deleteIfExists(output);}
		}
	}

	/** Wraps a fresh dump; called only after native validation by pack(). */
	private static void envelope(final String dump, final String output, final byte[][] provenance) throws IOException{
		final ByteFile in=ByteFile.makeByteFile(dump, false);
		final ByteStreamWriter out=new ByteStreamWriter(output, true, false, false);
		out.start();
		final MessageDigest root=digest();
		MessageDigest family=null;
		boolean ended=false;
		try{
			final byte[] first=in.nextLine();
			if(first==null || !new String(first, StandardCharsets.US_ASCII).equals("#format\thbm_text_v1")){
				throw new IOException("Expected native dense dump header");
			}
			emit(out, root, "#format\t"+FORMAT);
			emit(out, root, "#graph_contract\t"+CONTRACT);
			for(int i=0; i<provenance.length; i++){emit(out, root, "#provenance_"+i+"\t"+DigestSuffix.fromDigest(provenance[i]));}
			for(byte[] line=in.nextLine(); line!=null; line=in.nextLine()){
				if(ended || line.length==0){throw new IOException("Empty/trailing dump row");}
				if(line[0]=='#'){
					if(family!=null){throw new IOException("Header inside family");}
					update(root, line); out.println(line);
				}else if(line.length==1 && line[0]=='e'){
					if(family==null){throw new IOException("Unmatched family end");}
					emit(out, root, "e\t"+DigestSuffix.fromDigest(family.digest())); family=null;
				}else if(line.length==1 && line[0]=='z'){
					if(family!=null){throw new IOException("Unterminated family");}
					emit(out, null, "z\t"+DigestSuffix.fromDigest(root.digest())); ended=true;
				}else{
					if(line[0]=='f'){
						if(family!=null){throw new IOException("Nested family");}
						family=digest();
					}
					if(family==null){throw new IOException("Body outside family");}
					update(family, line); out.println(line);
				}
			}
			if(!ended){throw new IOException("Missing dump terminator");}
		}finally{
			final boolean readError=in.close(), writeError=out.poisonAndWait();
			if(readError || writeError){throw new IOException("Dense envelope I/O failure");}
		}
	}

	private static void emit(final ByteStreamWriter out, final MessageDigest digest, final String line){
		final byte[] bytes=line.getBytes(StandardCharsets.UTF_8);
		update(digest, bytes); out.println(bytes);
	}
	private HbmDenseTextPacker(){}
}
