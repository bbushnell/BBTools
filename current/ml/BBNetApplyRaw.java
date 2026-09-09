package ml;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import shared.Tools;

/**
 * Gate harness: loads a {@code .bbnet} with the real production parser
 * ({@link CellNetParser}) and applies it to RAW vector rows (whitespace-separated
 * floats, one vector per line, no {@code #dims} header), printing predictions to
 * stdout in the GPU trainer's {@code --dump-preds} format so the two files diff
 * directly: a {@code #row<TAB>pred0[<TAB>pred1...]} header, then one line per row,
 * {@code rowIndex<TAB>%.8f...}. This is the Java half of the checkpoint-export
 * round-trip gate ({@code bbnet_export.py}); {@link BBNetCompare} is the strict
 * comparator that consumes both files.
 *
 * <p>Packaged copy of bbtools-dev {@code magqc_gpu_sweep/BBNetApply.java} (default
 * package, SHA-256 {@code 90c10fdd...} at bbtools-dev {@code c26965b}), renamed so it
 * does not collide with {@link BBNetApply} (the {@code #dims}-TSV runner). Behavior on
 * valid input is identical byte for byte; the only changes are the package, the
 * usage failure (throws instead of {@code System.exit}), an unreadable-file preflight
 * (throws instead of the parser's hard exit), a reader {@code finally}, and a stdout
 * error check at the end.</p>
 *
 * <p>Usage: {@code java -ea ml.BBNetApplyRaw <net.bbnet> <raw_rows.tsv>}</p>
 *
 * @author Ady (bbtools-dev, 2026-09), UMP45 (packaging, 2026-09-08)
 */
public class BBNetApplyRaw {

	public static void main(String[] args) throws Exception{
		if(args.length!=2){
			throw new IllegalArgumentException("Usage: java -ea ml.BBNetApplyRaw <net.bbnet> <raw_rows.tsv>");
		}
		//Preflight: a missing net is a hard JVM exit inside the parser's file layer
		//(ReadWrite.java:1345 exceptionKill), not an exception; reject it as one first.
		if(!new File(args[0]).canRead()){throw new IllegalArgumentException("net file not readable: "+args[0]);}
		if(!new File(args[1]).canRead()){throw new IllegalArgumentException("rows file not readable: "+args[1]);}
		final CellNet net=CellNetParser.load(args[0], false);
		if(net==null){throw new RuntimeException("failed to load "+args[0]);}
		final BufferedReader br=new BufferedReader(new FileReader(args[1]));
		String header=null;
		int row=0;
		try{
			for(String line=br.readLine(); line!=null; line=br.readLine()){
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				final String[] parts=line.trim().split("\\s+");
				final float[] vec=new float[parts.length];
				for(int i=0; i<parts.length; i++){vec[i]=Float.parseFloat(parts[i]);}
				net.applyInput(vec);
				net.feedForward();
				final float[] out=net.getOutput();
				if(header==null){
					final StringBuilder hb=new StringBuilder("#row");
					for(int k=0; k<out.length; k++){hb.append('\t').append("pred").append(k);}
					System.out.println(hb);
					header=hb.toString();
				}
				final StringBuilder sb=new StringBuilder();
				sb.append(row);
				for(float v : out){sb.append('\t').append(Tools.format("%.8f", v));}
				System.out.println(sb);
				row++;
			}
		}finally{
			br.close();
		}
		if(row==0){throw new RuntimeException(args[1]+": zero data rows");}
		if(System.out.checkError()){throw new RuntimeException("BBNetApplyRaw: error writing predictions to stdout");}
	}
}
