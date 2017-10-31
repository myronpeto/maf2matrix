package org.panda.resource;

import org.panda.resource.tcga.MutTuple;
import org.panda.resource.tcga.MutationReader;
// import org.panda.utility.ArrayUtil;
//import org.panda.utility.StringUtil;
// import org.panda.utility.statistics.FDR;
// import org.panda.resource.tcga.AlterationMatrixSeparator;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Converts the PanCan MAF to mutation matrix.
 *
 * @author Ozgun Babur
 */
public class PanCanMAFToMatrix
{
	public static void main(String[] args) throws IOException
	{
		if (args.length != 2) {
			System.err.println("usage: PanCanMAFToMatrix maf-file out-dir");
			System.exit(1);
		}
		convertToMatrix(args[0], args[1]);
//		convertToMatrixWithSelectMutations();
//		separateToChunks();
	}


//	public static void separateToChunks() throws IOException
//	{
//		for (int k = 7; k <= 10; k++)
//		{
//			String base = "/home/babur/Documents/mutex/TCGA/PanCan-shuffled-" + k + "/";
//			System.out.println("base = " + base);
//			String wholeDir = base + "1/1";
//
//			for (int i = 2; i <= 30; i++)
//			{
//				String dir = base + i;
//				Files.createDirectories(Paths.get(dir));
//
//				AlterationMatrixSeparator.separate(dir, wholeDir, i);
//			}
//		}
//	}

	public static void convertToMatrix(String pancanMAF, String outDir) throws IOException
	{
		System.out.println("maf file: ");
		System.out.println(pancanMAF);

		MutationReader mr = new MutationReader(pancanMAF, "Missense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
			"Nonsense_Mutation", "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Translation_Start_Site");

		System.out.println("Getting samples");

		String[] samples = mr.getSamples().stream().sorted().collect(Collectors.toList()).toArray(new String[0]);


		System.out.println("Writing to output file");
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outDir + "/" + pancanMAF.substring(0, pancanMAF.length()-4) + "_Matrix.txt"));
		for (String sample : samples)
		{
			writer.write("\t" + sample);
		}
		for (String gene : mr.getGenes())
		{
			boolean[] mut = mr.getGeneAlterationArray(gene, samples);
			if (mut != null)
			{
				writer.write("\n" + gene);
				for (boolean m : mut)
				{
					writer.write("\t" + (m ? "1" : "0"));
				}
			}
		}

		System.out.println(pancanMAF.substring(0, pancanMAF.length()-4) + "_Matrix.txt written.");
		writer.close();
	}

//	private static boolean[] getAlterations(List<MutTuple>[] muts, Set<Integer> hotLocs, boolean isTrunBiased)
//	{
//		boolean[] b = new boolean[muts.length];
//
//		for (int i = 0; i < b.length; i++)
//		{
//			boolean m = false;
//			if (muts[i] != null) // we don't handle no-data cases
//			{
//				for (MutTuple mut : muts[i])
//				{
//					if (mut.isDeleterious() && isTrunBiased)
//					{
//						m = true;
//						break;
//					}
//
//					if (hotLocs != null)
//					{
//						int loc = StringUtil.findFirstInt(mut.value);
//						if (hotLocs.contains(loc))
//						{
//							m = true;
//							break;
//						}
//					}
//				}
//			}
//			b[i] = m;
//		}
//		return b;
//	}
//
}
