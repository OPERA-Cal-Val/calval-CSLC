#!/usr/bin/env python3

"""Automated routine to prepare a CSLC Cal/Val report"""
import argparse
from pathlib import Path

import PyPDF2
import pandas as pd
from PIL import Image as PIL_Image
import matplotlib.pyplot as plt
from matplotlib import rc
from tqdm import tqdm

from wand.color import Color
from wand.drawing import Drawing
from wand.image import Image

# set matplotlib font
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

def create_parser():
    """ Initiate input parser"""

    parser = argparse.ArgumentParser(description='Create CSLC '
                                     'validation report.')
    parser.add_argument('-p', '--pdir', dest='pdir', type=str,
                        help='Specify directory containing pngs')
    parser.add_argument('-c', '--cdir', dest='cdir', type=str,
                        help='Specify directory containing csvs')
    parser.add_argument('-o', '--outdir', dest='outdir', type=str,
                        default='./', help='Specify output directory')
    parser.add_argument('-f', '--fname', dest='fname', type=str,
                        default='output.pdf', help='Specify output filename')

    return parser


def cmd_line_parse(iargs=None):
    """ Wrap around input parser"""
    parser = create_parser()
    inputs = parser.parse_args(args=iargs)

    return inputs


def merge_pdf_pages(input_pdf,
                    output_pdf,
                    width,
                    height,
                    resolution,
                    top=0,
                    left=0):
    """ Merge two PDF pages into one"""
    # Create a new image (blank page)
    with Image(width=width,
               height=(height+top)*2,
               resolution=resolution) as new_page:
        # Read the first page of the input PDF
        with Image(filename=f'{input_pdf}[0]') as page1:
            # Set the dimensions for the new page based on the first page
            new_page.width = page1.width
            new_page.height = page1.height * 2  # Concat 2 pages vertically

            # Composite the first page onto the new page
            new_page.composite(page1, top=top, left=left)

            # Read the second page of the input PDF
            with Image(filename=f'{input_pdf}[1]') as page2:
                # Composite the second page below the first page
                new_page.composite(page2, top=page1.height+top, left=0)

        # add border (color, width, and height)
        new_page.border('White', 12, 0)
        # Save the concatenated PDF
        new_page.save(filename=output_pdf)

    return


def concatenate_pdf_pages(input_pdf,
                          output_pdf):
    """ Consolidate two PDFs into one"""
    # Make copy of input
    input_pdf1 = str(output_pdf).replace('.pdf', '_temp.pdf')
    Path(output_pdf).replace(input_pdf1)

    # Open the input PDF files
    with open(input_pdf1, 'rb') as file1, open(input_pdf, 'rb') as file2:
        # Create PDF reader objects
        pdf_reader1 = PyPDF2.PdfFileReader(file1)
        pdf_reader2 = PyPDF2.PdfFileReader(file2)

        # Create a PDF writer object for the output file
        pdf_writer = PyPDF2.PdfFileWriter()

        # Add all pages from the first PDF
        for page_num in range(pdf_reader1.numPages):
            page = pdf_reader1.getPage(page_num)
            pdf_writer.addPage(page)

        # Add all pages from the second PDF
        for page_num in range(pdf_reader2.numPages):
            page = pdf_reader2.getPage(page_num)
            pdf_writer.addPage(page)

        # Save the concatenated PDF to the output file
        with open(output_pdf, 'wb') as output_file:
            pdf_writer.write(output_file)

    # delete temp file
    Path(input_pdf1).unlink()

    return


def add_text(output_pdf,
             header_txt,
             font_size=48,
             fill_color='black',
             pos_x=40,
             pos_y=60):
    """ Add title text for a given page"""
    # re-open cal/val report to add text to top of page
    pdf_tot = Image(filename=output_pdf)

    # Create a Drawing content to add to image
    with Drawing() as draw:
        # Set the font properties
        draw.font = 'Helvetica'
        draw.font_size = font_size
        draw.fill_color = Color(fill_color)
        # add gray background to text
        draw.text_under_color = Color('#CCCCCC')
        # Add the text to the image
        draw.text(pos_x, pos_y, header_txt)
        draw(pdf_tot)

    # Save the image with the added text
    pdf_tot.save(filename=output_pdf)

    return


def resize_img(source_file,
               width,
               height,
               resolution):
    """ Resize an input image or PDF page"""
    # Open image and set dimensions/res
    outer_img = Image(width=width,
                      height=height,
                      resolution=resolution,
                      units='pixelsperinch')

    # if larger than len x wid, fit within box, preserving aspect ratio
    # accordingly determine output dimensions
    tmp_img = Image(filename=source_file)
    tmp_img.transform(resize=f'{width}x{height}>')
    src_wid = int(tmp_img.size[0] * .92)
    src_hgt = int(tmp_img.size[1] * .92)

    del tmp_img

    # Use different libraries to resample PDFs and PNGs
    tmp_source_file = str(source_file)[:-3] + 'temp.png'
    if str(source_file)[-3:] == 'pdf':
        img = Image(filename=source_file)
    else:
        img = PIL_Image.open(source_file)
        # Resize the image using the Lanczos filter (high-quality)
        resized_img = img.resize((src_wid, src_hgt), PIL_Image.LANCZOS)
        # Save the resized image
        resized_img.save(tmp_source_file, dpi=(1200, 1200))
        del img, resized_img
        img = Image(filename=tmp_source_file)
        Path(tmp_source_file).unlink()
    outer_img.composite(img,
                        left=int((width - img.width) / 2),
                        top=int((height - img.height) / 2))

    return outer_img


def png_wrapper(pngs_paired,
                burst,
                output_pdf,
                width,
                height,
                resolution,
                visual_type,
                top_buffer=0,
                fill_color='black'):
    """ Loop through all input images"""
    # loop through each summary png
    for input_pngs in pngs_paired:
        if burst in input_pngs[0]:
            # build temporary file
            temp_pdf = input_pngs[0].replace('.png', '_temp.pdf')
            for png in input_pngs:
                update_pdf(temp_pdf, png, width, height, resolution)

            if len(input_pngs) > 1:
                # combine images between 2 pages into one
                temp_pdf2 = input_pngs[0].replace('.png', '_fin.pdf')
                merge_pdf_pages(temp_pdf, temp_pdf2,
                                width,
                                height,
                                resolution,
                                top=top_buffer)
            else:
                temp_pdf2 = input_pngs[0].replace('.png', '_temp.pdf')

            # rescale pages
            pdf_rescale = resize_img(temp_pdf2,
                                     width,
                                     height*2,
                                     resolution)
            pdf_rescale.save(filename=temp_pdf2)
            del pdf_rescale

            # add text
            if visual_type == 'ALE_summary':
                header_txt = f'ALE summary for burst {burst}'
            if visual_type == 'ALE_indiv':
                burst_date = input_pngs[0]
                burst_date = burst_date.split('.png')[0]
                burst_date = burst_date.split('_')[-1]
                header_txt = f'ALE: burst {burst}, date {burst_date}'
            if visual_type == 'RLE_indiv':
                header_txt = f'RLE: burst {burst}'
            add_text(temp_pdf2, header_txt, font_size=32,
                     fill_color=fill_color, pos_x=40, pos_y=40)

            # append page to cal/val report
            concatenate_pdf_pages(temp_pdf2, output_pdf)

            # delete temp files
            if Path(temp_pdf).exists():
                Path(temp_pdf).unlink()
            if Path(temp_pdf2).exists():
                Path(temp_pdf2).unlink()

    return


def csv_wrapper(input_csv,
                output_pdf):
    """ Loop through all input csvs"""
    # Read the CSV data using Pandas
    df = pd.read_csv(input_csv, dtype=str)

    # assign 3 sig figs to measurements
    df_val_cols = ['Delta Ground Range (m)',
                   'Delta Ground Range stdev (m)',
                   'Delta Azimuth (m)',
                   'Delta Azimuth stdev (m)']
    for i in df_val_cols:
        df[i] = \
            df[i].astype(float).apply(lambda x:
                                      round(x, 3)).astype(str)

    # Combine stdev with +/- symbol
    df['\u0394 Ground Range (m)'] = \
        df['Delta Ground Range (m)'].astype(str) + \
        '\u00B1' + \
        df['Delta Ground Range stdev (m)'].astype(str)
    df.drop(['Delta Ground Range (m)', 'Delta Ground Range stdev (m)'],
            axis=1, inplace=True)
    df['\u0394 Azimuth (m)'] = \
        df['Delta Azimuth (m)'].astype(str) + \
        '\u00B1' + \
        df['Delta Azimuth stdev (m)'].astype(str)
    df.drop(['Delta Azimuth (m)', 'Delta Azimuth stdev (m)'],
            axis=1, inplace=True)

    # highlight failed row(s) in pink
    table_wid = df.shape[1]
    if 'ALE' in Path(output_pdf).name:
        cell_colors = [i for i in df['ALE (Pass/Fail)']]
    if 'RLE' in Path(output_pdf).name:
        cell_colors = [i for i in df['RLE (Pass/Fail)']]
    for i in enumerate(cell_colors):
        if i[1].lower() == 'pass':
            cell_colors[i[0]] = ['white'] * table_wid
        else:
            cell_colors[i[0]] = ['#FFC0CB'] * table_wid

    # set dims, and scale by expected page size (hardcoded ratios)
    fig_w = 16 * (800/576)
    fig_l = 24 * (600/432)

    # Create a figure and axis for the plot
    fig = plt.figure(figsize=(fig_w, fig_l))
    ax = fig.add_subplot(111)

    # Create a table from the Pandas DataFrame
    table = ax.table(cellText=df.values,
                     colLabels=df.columns,
                     loc='center',
                     cellColours=cell_colors)

    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(30)
    table.scale(1.2, 2.8)

    # Remove the axis
    ax.axis('off')

    # Save the figure as a PDF
    plt.savefig(output_pdf, pad_inches=0.1)

    # Close the figure
    plt.close()

    return


def update_pdf(output_pdf,
               png,
               width,
               height,
               resolution,
               resize=True):
    """ Update a PDF with new content"""
    # initialize output file
    if resize is True:
        pdf_init = resize_img(png,
                              width,
                              height,
                              resolution)
    else:
        pdf_init = Image(filename=png)

    # Save the PDF
    # Initiate output if it does not exist
    if not Path(output_pdf).exists():
        pdf_init.save(filename=output_pdf)
    # otherwise, create a temp file
    else:
        temp_pdf = png.replace('.png', '.pdf')
        # create temp file
        if resize is True:
            pdf_init.save(filename=temp_pdf)
        # re-open cal/val report to append each new page
        concatenate_pdf_pages(temp_pdf, output_pdf)
        # delete temp file
        if resize is True:
            Path(temp_pdf).unlink()

    del pdf_init

    return


def prep_cslc_report(inps):
    """ Create a CSLC Cal/Val report"""
    # hardcode output page dims/resolution
    width = 1600
    height = 1200
    resolution = 1200

    # pass inputs as Path objects
    inps.pdir = Path(inps.pdir)
    inps.cdir = Path(inps.cdir)
    inps.outdir = Path(inps.outdir)

    # if necessary create output directory
    if not inps.outdir.exists():
        inps.outdir.mkdir(parents=True, exist_ok=True)

    # check if specified output file already exists
    output_pdf = inps.outdir / inps.fname
    if output_pdf.exists():
        raise Exception(f'Specified output -f {inps.fname} already exists')

    # find all pngs
    png_files = [str(i) for i in inps.pdir.glob('*.png')]
    if not png_files:
        raise Exception(f'No images found in specified path -p {inps.pdir}')

    # exit if a mix of ALE and RLE figures are found in a given subdirectory
    ale_pngs = [i for i in png_files if 'ALE' in Path(i).name]
    rle_pngs = [i for i in png_files if 'RLE' in Path(i).name]
    if ale_pngs != [] and rle_pngs != []:
        raise Exception(f'Mix of RLE and ALE pngs found in -p {inps.pdir}')

    # separate out summary pngs
    summary_pngs = [i for i in png_files if 'summary' in Path(i).name]
    ale_pdf = False
    # sort and split into pairs, if ALE figures captured
    if summary_pngs != []:
        ale_pdf = True
        summary_pngs = sorted(summary_pngs,
                              key=lambda x: Path(x).name[:38])
        summary_pngs_paired = [summary_pngs[i:i + 2]
                               for i in range(0, len(summary_pngs), 2)]

    # get individual burst pngs
    burst_pngs = [i for i in png_files if i not in summary_pngs]
    # sort and split into pairs, if ALE figures captured
    if summary_pngs != []:
        burst_pngs = sorted(burst_pngs, key=lambda x:
                            Path(x).name[:3], reverse=True)
        burst_pngs = sorted(burst_pngs,
                            key=lambda x: x[-28:])
        burst_pngs_paired = [burst_pngs[i:i + 2]
                             for i in range(0, len(burst_pngs), 2)]
        # get list of burst IDs
        burst_ids = [Path(i[0]).name[26:-13]
                     for i in burst_pngs_paired]
    # sort RLE figures
    else:
        burst_pngs_paired = sorted(burst_pngs,
                                   key=lambda x: x[-19:])
        burst_pngs_paired = [[i] for i in burst_pngs_paired]
        burst_ids = [Path(i[0]).name.split('.png')[0][4:]
                     for i in burst_pngs_paired]

    # get list of burst IDs
    burst_ids = sorted(list(set(burst_ids)))


    # find all csvs
    csvs_files = [str(i) for i in inps.cdir.glob('*.csv')]
    if not csvs_files:
        raise Exception(f'No csvs found in specified path -c {inps.cdir}')

    # exit if a mix of ALE and RLE csvs are found in a given subdirectory
    ale_csvs = [i for i in csvs_files if 'ALE' in Path(i).name]
    rle_csvs = [i for i in csvs_files if 'RLE' in Path(i).name]
    if ale_csvs != [] and rle_csvs != []:
        raise Exception(f'Mix of RLE and ALE csvs found in -p {inps.cdir}')


    # initialize output file with validation map
    validation_map = 'cslc_validation_sites.png'
    update_pdf(output_pdf, validation_map, width, height, resolution)
    # add text
    header_txt = 'Validation sites'
    add_text(output_pdf, header_txt, pos_x=40, pos_y=40)
    # rescale page
    pdf_rescale = resize_img(output_pdf,
                             width,
                             height*2,
                             resolution)
    pdf_rescale.save(filename=output_pdf)
    del pdf_rescale


    # add table from CSV file
    for csv_f in csvs_files:
        output_csvpdf = Path(csv_f).name.replace('.csv', '.pdf')
        output_csvpdf = inps.outdir / output_csvpdf
        # create pdf from table
        csv_wrapper(csv_f, output_csvpdf)
        # append to report
        concatenate_pdf_pages(output_csvpdf, output_pdf)
        # delete temp file
        output_csvpdf.unlink()


    # loop through each burst
    for burst in tqdm(burst_ids, desc="Processing burst"):
        # if ALE figures captured
        print(burst)
        if ale_pdf is True:
            # loop through each summary png
            png_wrapper(summary_pngs_paired, burst, output_pdf,
                        width, height, resolution,
                        visual_type='ALE_summary',
                        top_buffer=60, fill_color='red')
            # loop through each ALE png
            png_wrapper(burst_pngs_paired, burst, output_pdf,
                        width, height, resolution,
                        visual_type='ALE_indiv',
                        top_buffer=0, fill_color='black')
        else:
            # loop through each ALE png
            png_wrapper(burst_pngs_paired, burst, output_pdf,
                        width, height, resolution,
                        visual_type='RLE_indiv',
                        top_buffer=0, fill_color='black')

    print(f'PDF created: {output_pdf}')

    return


if __name__ == '__main__':
    calval_inps = cmd_line_parse()
    prep_cslc_report(calval_inps)
