const http = require('http');
const fs = require('fs');
const path = require('path');

const PORT = 8090;
const ROOT = path.join(__dirname, '..');

const MIME = {
    '.html': 'text/html',
    '.js': 'application/javascript',
    '.css': 'text/css',
    '.json': 'application/json',
    '.png': 'image/png',
    '.jpg': 'image/jpeg',
    '.jfem': 'application/octet-stream',
    '.md': 'text/plain',
};

http.createServer((req, res) => {
    let filePath = path.join(ROOT, req.url === '/' ? 'POST/postv11.html' : req.url);
    const ext = path.extname(filePath).toLowerCase();
    fs.readFile(filePath, (err, data) => {
        if (err) { res.writeHead(404); res.end('Not found'); return; }
        res.writeHead(200, { 'Content-Type': MIME[ext] || 'application/octet-stream' });
        res.end(data);
    });
}).listen(PORT, () => console.log(`Serving on http://localhost:${PORT}`));
