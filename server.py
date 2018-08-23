#!/usr/bin/python

from flask import Flask, render_template, request
import scipy.stats

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/results', methods=['GET', 'POST'])
def results():
    foreground_organisms = set()
    if 'foreground' in request.form:
        for organism in request.form['foreground'].split('\n'):
            organism = organism.strip()
            if organism:
                foreground_organisms.add(organism)
    background_organisms = set()
    if 'background' in request.form:
        for organism in request.form['background'].split('\n'):
            organism = organism.strip()
            if organism:
                background_organisms.add(organism)
    dictionary = 'disease'
    if 'dictionary' in request.form:
        dictionary = request.form['dictionary']
    z_cutoff = 5.0
    if 'zcutoff' in request.form:
        z_cutoff = float(request.form['zcutoff'])
    p_cutoff = 1.0
    if 'pcutoff' in request.form:
        p_cutoff = float(request.form['pcutoff'])

    # Read and count associations above cutoff for foreground and background.
    foreground_counts = {}
    background_counts = {}
    term_names = {}
    for line in open('/data/download/organism_%s_textmining_full.tsv' % dictionary, 'r'):
        fields = line.strip().split('\t')
        z_score = fields[4]
        if float(z_score) >= z_cutoff:
            organism = fields[0]
            term = fields[2]
            name = fields[3]
            if organism in foreground_organisms:
                if term in foreground_counts:
                    foreground_counts[term] += 1
                else:
                    foreground_counts[term] = 1
                    if term not in term_names:
                        term_names[term] = name
                if organism in background_organisms:
                    if term in background_counts:
                        background_counts[term] += 1
                    else:
                        background_counts[term] = 1

    # Loop over counted terms, calculate statistics, and print results.
    results = []
    for term in foreground_counts.keys():
        name = term_names[term]
        foreground_count = foreground_counts[term]
        if background_organisms:
            background_count = 0
            if term in background_counts:
                background_count = background_counts[term]-foreground_count
                foreground_size = len(foreground_organisms)
                background_size = len(background_organisms)-foreground_size
                table = [[foreground_count,foreground_size-foreground_count],[background_count,background_size-background_count]]
                ratio, p_value = scipy.stats.fisher_exact(table, alternative='greater')
                if p_value < p_cutoff:
                    results.append((term, name, foreground_count, background_count, '%.9f' % p_value))
        else:
            results.append((term, name, foreground_count))

    return render_template('results.html', results=results, enrichment=bool(background_counts))

if __name__ == '__main__':
    app.run(debug=False, host='0.0.0.0')
